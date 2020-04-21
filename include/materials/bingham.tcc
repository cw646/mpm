//! Constructor with material properties
template <unsigned Tdim>
mpm::Bingham<Tdim>::Bingham(unsigned id, const Json& material_properties)
    : Material<Tdim>(id, material_properties) {
  try {
    density_ = material_properties.at("density").template get<double>();
    youngs_modulus_ =
        material_properties.at("youngs_modulus").template get<double>();
    poisson_ratio_ =
        material_properties.at("poisson_ratio").template get<double>();
    tau0_ = material_properties.at("tau0").template get<double>();
    mu_ = material_properties.at("mu").template get<double>();
    critical_shear_rate_ =
        material_properties["critical_shear_rate"].template get<double>();
    // Calculate bulk modulus
    bulk_modulus_ = youngs_modulus_ / (3.0 * (1. - 2. * poisson_ratio_));

    properties_ = material_properties;
  } catch (Json::exception& except) {
    console_->error("Material parameter not set: {} {}\n", except.what(),
                    except.id);
  }
}

//! Initialise history variables
//CW comment: some state vars will eventually be added via input file.
template <unsigned Tdim>
mpm::dense_map mpm::Bingham<Tdim>::initialise_state_variables() {
  mpm::dense_map state_vars = {{"pressure", 0.0},
                               {"floc", 0.0},
                               {"rest_t", 0.0},
                               {"alpha", 0.0116},
                               {"a_thix", 0.7},
                               {"dt",1E-4},
			       {"shear_rate", 0.0},
			       {"visco", 0.0},
			       {"strain0", 0.0},
			       {"strain1", 0.0},
			       {"strain2", 0.0},
			       {"strain3", 0.0},
			       {"strain4", 0.0},
			       {"strain5", 0.0}};
  return state_vars;
}

//! Compute pressure
template <unsigned Tdim>
double mpm::Bingham<Tdim>::thermodynamic_pressure(
    double volumetric_strain) const {
  // Bulk modulus
  return (-bulk_modulus_ * volumetric_strain);
}

//! Compute stress
template <unsigned Tdim>
Eigen::Matrix<double, 6, 1> mpm::Bingham<Tdim>::compute_stress(
    const Vector6d& stress, const Vector6d& dstrain,
    const ParticleBase<Tdim>* ptr, mpm::dense_map* state_vars) {


  // Get strain rate
  auto strain_rate = ptr->strain_rate();

  // Convert strain rate to rate of deformation tensor
  strain_rate.tail(3) *= 0.5;

  // Set threshold for minimum critical shear rate
  const double shear_rate_threshold = 1.0E-15;
  if (critical_shear_rate_ < shear_rate_threshold)
    critical_shear_rate_ = shear_rate_threshold;

  // Rate of shear = sqrt(2 * D_ij * D_ij)
  // Since D (D_ij) is in Voigt notation (D_i), and the definition above is in
  // matrix, the last 3 components have to be doubled D_ij * D_ij = D_0^2 +
  // D_1^2 + D_2^2 + 2*D_3^2 + 2*D_4^2 + 2*D_5^2 Yielding is defined: rate of
  // shear > critical_shear_rate_^2 Checking yielding from strain rate vs
  // critical yielding shear rate
  double shear_rate =
      std::sqrt(2. * (strain_rate.dot(strain_rate) +
                      strain_rate.tail(3).dot(strain_rate.tail(3))));

// Apparent_viscosity maps shear rate to shear stress

  double apparent_viscosity = 0.;

  Eigen::Matrix<double, 6, 1> tau;

  double tollerance = 0.0001;

  (*state_vars)["shear_rate"] = shear_rate;

  (*state_vars)["strain0"] = strain_rate(0);
  (*state_vars)["strain1"] = strain_rate(1);
  (*state_vars)["strain2"] = strain_rate(2);
  (*state_vars)["strain3"] = strain_rate(3);
  (*state_vars)["strain4"] = strain_rate(4);
  (*state_vars)["strain5"] = strain_rate(5);

//console_->info("Strain1 {}  ", (strain_rate(2,2)));

  if (shear_rate * shear_rate > critical_shear_rate_ * critical_shear_rate_) {

    if ((*state_vars)["floc"] > tollerance) {  // is lambda greater than 0

      //console_->info("shear_rate {} ", shear_rate);

      auto floc_prev = (*state_vars)["floc"];  // previous floculation state

      double dldt = -(*state_vars)["alpha"] * (*state_vars)["floc"] *
                    shear_rate;             // calculate change in lambda

      auto floccheck = (*state_vars)["floc"] + dldt;  // check to see if lambda positive

      if (floccheck > tollerance) {      // is new lambda > 0

        (*state_vars)["floc"] = floc_prev + dldt;  // update floculation state

        // current floculation state
        double tau0_temp =
            (1 + (*state_vars)["floc"]) * tau0_;  // caluclate new apparent yield stress


        apparent_viscosity = 2. * ((tau0_temp / shear_rate) + mu_);


        // apparent rest time
        (*state_vars)["rest_t"] = (tau0_temp - tau0_) / (*state_vars)["a_thix"];

      } else {
        (*state_vars)["floc"] = 0.;

        apparent_viscosity = 2. * ((tau0_ / shear_rate) + mu_);

        (*state_vars)["rest_t"] = 0;
      }

    } else {

      apparent_viscosity = 2. * ((tau0_ / shear_rate) + mu_);
      (*state_vars)["rest_t"] = 0;

    }


    tau = apparent_viscosity * strain_rate;

  } else {

    apparent_viscosity = 1000*mu_;

    tau = (2 * apparent_viscosity) * strain_rate;
    (*state_vars)["rest_t"] += (*state_vars)["dt"];

    //(*state_vars)["tau0"] = tau_t;

    double tau_t = tau0_ + ((*state_vars)["a_thix"] * (*state_vars)["rest_t"]);

    (*state_vars)["floc"] = (tau_t / tau0_) - 1;  // update floculation state with
                                        // increase in thixotropic presence

    //console_->info("floc {} ", (*state_vars)["floc"]);
  }

  (*state_vars)["visco"] = apparent_viscosity;

  // Update pressure
  (*state_vars).at("pressure") +=
      (compressibility_multiplier_ *
       this->thermodynamic_pressure(ptr->dvolumetric_strain()));

  // Update volumetric and deviatoric stress
  // thermodynamic pressure is from material point
  // stress = -thermodynamic_pressure I + tau, where I is identity matrix or
  // direc_delta in Voigt notation
  const Eigen::Matrix<double, 6, 1> updated_stress =
      -(*state_vars).at("pressure") * this->dirac_delta() *
          compressibility_multiplier_ +
      tau;

  return updated_stress;
}

//! Dirac delta 2D
template <>
inline Eigen::Matrix<double, 6, 1> mpm::Bingham<2>::dirac_delta() const {

  return (Eigen::Matrix<double, 6, 1>() << 1.f, 1.f, 0.f, 0.f, 0.f, 0.f)
      .finished();
}

//! Dirac delta 3D
template <>
inline Eigen::Matrix<double, 6, 1> mpm::Bingham<3>::dirac_delta() const {

  return (Eigen::Matrix<double, 6, 1>() << 1.f, 1.f, 1.f, 0.f, 0.f, 0.f)
      .finished();
}
