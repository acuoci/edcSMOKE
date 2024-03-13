/*-----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Author: Alberto Cuoci <alberto.cuoci@polimi.it>                       |
|   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of OpenSMOKE++ framework.                           |
|                                                                         |
|   License                                                               |
|                                                                         |
|   Copyright(C) 2021  Alberto Cuoci                                      |
|   OpenSMOKE++ is free software: you can redistribute it and/or modify   |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   OpenSMOKE++ is distributed in the hope that it will be useful,        |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with OpenSMOKE++. If not, see <http://www.gnu.org/licenses/>.   |
|                                                                         |
\*-----------------------------------------------------------------------*/

#ifndef OpenSMOKE_LookUpTable_ZY_H
#define OpenSMOKE_LookUpTable_ZY_H

// External parameters
#include <Eigen/Dense>

// OpenSMOKE++ Definitions
#include "OpenSMOKEpp"

// CHEMKIN maps
#include "maps/Maps_CHEMKIN"


namespace OpenSMOKE
{
	//!  A class to manage lookup tables based on the mixture-fraction/progress-variable approach
	/*!
	A class to manage lookup tables based on the mixture-fraction/progress-variable approach
	*/

	class LookUpTable_ZY
	{

	public:

		enum Initialization_Types { FROM_SCRATCH, FROM_BACKUP, FROM_RECONSTRUCTION, FROM_ONTHEFLY_SOLUTION };

	private:

		enum LookUpTable_ZY_Type { UNIFORM_Z_UNIFORM_Y };
		enum Interpolated_Types { TEMPERATURE, DENSITY, SPECIFIC_HEAT, THERMAL_DIFFUSIVITY, MASS_FRACTION };

	public:

		/**
		*@brief Default constructor
		*@param thermodynamicsMap	reference to the thermodynamic map
		*/
		LookUpTable_ZY(OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap);

		/**
		*@brief Setup from a dictionary
		*@param dictionary dictionary name
		*/
		void SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary);

		/**
		*@brief Setup from user provided path
		*@param path_to_xml path to the XML lookup table
		*@param is_apriori if true, a priori analysis is performed
		*/
		void Setup(const boost::filesystem::path& path_to_xml, const std::vector<std::string>& list_fields, const bool is_apriori);

		/**
		*@brief Interpolate
		*@param x mixture-fraction (-)
		*@param c progress-variable (-)
		*/
		void Interpolate(const double z, const double c);

		/**
		*@brief Interpolate a specific variable, assuming the production/consumption terms are available
		*@param k index of variable
		*@param x mixture-fraction (-)
		*@param c progress-variable (-)
		*/
		void Interpolate(const unsigned int k, const double z, const double c);

		/**
		*@brief Reconstruct the mixture fraction from the mass fractions
		*@param Y mass fraction vector
		*@return the mixture fraction
		*/
		double ReconstructMixtureFraction(const double* Y);

		/**
		*@brief Reconstruct the progress variable from the mass fractions
		*@param Y mass fraction vector
		*@return the progress variable
		*/
		double ReconstructProgressVariable(const double* Y);

		/**
		*@brief Returns the interpolated mass fraction (according to the name provided as input argument)
		*@param name name of the required field
		*@return the interpolated mass fraction
		*/
		double Y(const std::string name) const;

		/**
		*@brief Returns the interpolated production term (according to the name provided as input argument)
		*@param name name of the required field
		*@return the interpolated production term (kg/m3/s)
		*/
		double OmegaP(const std::string name) const;

		/**
		*@brief Returns the interpolated consumption term (according to the name provided as input argument)
		*@param name name of the required field
		*@return the interpolated consumption term (kg/m3/s)
		*/
		double OmegaD(const std::string name) const;

		/**
		*@brief Returns the j-th interpolated mass fraction (according to the order reported in the list_fields variable)
		*@param j index of required interpolated value
		*@return the j-th interpolated mass fraction
		*/
		double Y(const unsigned int j) const { return interpolated_Y_[j]; }

		/**
		*@brief Returns the j-th interpolated production term (according to the order reported in the list_fields variable)
		*@param j index of required interpolated value
		*@return the j-th interpolated production term (kg/m3/s)
		*/
		double OmegaP(const unsigned int j) const { return interpolated_omegap_[j]; }

		/**
		*@brief Returns the j-th interpolated consumption term (according to the order reported in the list_fields variable)
		*@param j index of required interpolated value
		*@return the j-th interpolated consumption term (kg/m3/s)
		*/
		double OmegaD(const unsigned int j) const { return interpolated_omegad_[j]; }

		/**
		*@brief Returns the list of fields available in the lookup table
		*@return list of fields available in the lookup table
		*/
		const std::vector<std::string>& list_fields() const { return list_fields_; }

		/**
		*@brief Returns the indices of species belonging to the PAH12 family
		*@return indices of species belonging to the PAH12 family
		*/
		const std::vector<unsigned int>& indices_pah12() const { return indices_pah12_; }

		/**
		*@brief Returns the indices of species belonging to the PAH34 family
		*@return indices of species belonging to the PAH34 family
		*/
		const std::vector<unsigned int>& indices_pah34() const { return indices_pah34_; }

		/**
		*@brief Returns the indices of species belonging to the PAHLP family
		*@return indices of species belonging to the PAHLP family
		*/
		const std::vector<unsigned int>& indices_pahlp() const { return indices_pahlp_; }

		/**
		*@brief Returns the indices of species belonging to the SP family
		*@return indices of species belonging to the SP family
		*/
		const std::vector<unsigned int>& indices_sp() const { return indices_sp_; }

		/**
		*@brief Returns the indices of species belonging to the AGG family
		*@return indices of species belonging to the AGG family
		*/
		const std::vector<unsigned int>& indices_agg() const { return indices_agg_; }

		/**
		*@brief Returns true if the production/consumption terms are available for the given field
		*@return true if the production/consumption terms are available
		*/
		const std::vector<bool>& available_omega() const { return available_omega_; }

		/**
		*@brief Returns the indices of fields for which the production/consumption terms are available
		*@return indices of fields for which the production/consumption terms are available
		*/
		const std::vector<unsigned int>& list_available_omega() const { return list_available_omega_; }

		/**
		*@brief Returns true if the production/consumption terms are available for the given field
		*@param name the name of the field
		*@return true if the production/consumption terms are available
		*/
		bool available_omega(const std::string name) const;

		/**
		*@brief Returns the index of the requested species/class
		*@param name the requested species/class
		*@return index of the requested species/class (0-based)
		*/
		int Index(const std::string name) const;

		/**
		*@brief Returns the number of fields with production/consumption terms, i.e. the ones for which transport equations have ti be solved
		*@return the number of fields with production/consumption terms
		*/
		unsigned int n_eqs() const { return list_available_omega_.size(); }

		/**
		*@brief Returns true if the PAH consumption has to be accounted for
		*@return true if the PAH consumption has to be accounted for
		*/
		bool PAHConsumption() const { return is_pah_consumption_; }

		/**
		*@brief Returns the thermophoretic model to be applied
		*@return the thermophoretic model to be applied
		*/
		unsigned int thermophoretic_model() const { return thermophoretic_model_;  }

		/**
		*@brief Returns the Lewis number for fields to be transported
		*@return the Lewis number
		*/
		double Le() const { return Le_;  }

		/**
		*@brief Returns the initialization type
		*@return the initialization type
		*/
		Initialization_Types initialization() const { return initialization_; }


	private:

		/**
		*@brief Read lookup table from XML file
		*/
		void ReadTableFromXMLFile();

		/**
		*@brief Check input options
		*/
		void CheckInput();

		/**
		*@brief Summary on the screen
		*/
		void Summary();

		/**
		*@brief Returns the minimum local progress variable value
		*@param mixture fraction value
		*@return minimum local progress variable value
		*/
		double MinimumLocalProgressVariable(const double z);

		/**
		*@brief Returns the maximum local progress variable value
		*@param mixture fraction value
		*@return maximum local progress variable value
		*/
		double MaximumLocalProgressVariable(const double z);

		/**
		*@brief Interpolate from normalized progress variable
		*@param z mixture-fraction (-)
		*@param ctilde normalized progress-variable (-)
		*/
		void InterpolateFromNormalizedProgressVariable(const double z, const double ctilde);

	private:

		bool is_active_;						//!< true only if the table has been activated
		bool is_apriori_;						//!< a priori analysis
		bool is_pah_consumption_;				//!< is PAH consumption accounted for?
		unsigned int thermophoretic_model_;		//!< the thermophoretic model
		double Le_;								//!< the Lewis number for species to be transported

		OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap_;		//!< thermodynamics map

		LookUpTable_ZY_Type type_;		//!< type of dynamic boundary condition

		boost::filesystem::path path_output_;	//!< path to output folder

		Initialization_Types initialization_;		//!< initialization type

		unsigned int nz_;					//<! number of points along mixture fraction field
		unsigned int nc_;				//<! number of points along progress variable field

		std::vector<double> z_;				//<! mixture-fraction coordinates
		std::vector<double> c_;				//<! progress-variable coordinates
		std::vector<double> ctilde_;			//<! normalized progress-variable coordinates
		std::vector<double> c_min_;			//<! minimum progress variable at each mixture fraction
		std::vector<double> c_max_;			//<! maximum progress variable at each mixture fraction

		double min_z_;				//!< minimum mixture-fraction
		double max_z_;				//!< maximum mixture-fraction
		double min_overall_c_;			//!< minimum progress-variable (overall)
		double max_overall_c_;			//!< maximum progress-variable (overall)

		double dz_;					//<! mixture-fraction spacing (only in case of uniform grid)
		double dc_;					//<! progress-variable spacing (only in case of uniform grid)
		double dctilde_;				//<! normalized progress-variable spacing (only in case of uniform grid)

		std::vector<std::string> list_fields_;						//<! list of fields 
		std::vector< std::vector< std::vector<double> > > tables_Y_;			//<! lookup tables
		std::vector< std::vector< std::vector<double> > > tables_omegap_;		//<! lookup tables
		std::vector< std::vector< std::vector<double> > > tables_omegad_;		//<! lookup tables

		double query_z_;		//!< current mixture fraction for which the interpolation was carried out
		double query_c_;		//!< current progress-variable for which the interpolation was carried out

		std::vector<unsigned int> Y_def_indices_;	//!< indices of species defining the progress variable
		std::vector<double> Y_def_weights_;		//!< weights of species defining the progress variable (already scaled by the molecular weights)

		unsigned int jC_;		//!< index of C atom in the matrix containing atomic composition of species
		unsigned int jO_;		//!< index of O atom in the matrix containing atomic composition of species
		unsigned int jH_;		//!< index of H atom in the matrix containing atomic composition of species

		double WC_;		//!< molecular weight of C atom (in kg/kmol)
		double WO_;		//!< molecular weight of O atom (in kg/kmol)
		double WH_;		//!< molecular weight of H atom (in kg/kmol)

		Eigen::VectorXd Y_fuel_;	//!< mass fractions of fuel stream
		Eigen::VectorXd Y_ox_;		//!< mass fractions of oxidizer stream

		double ZstarFuel_;		//!< un-normalized mixture fraction on the fuel side
		double ZstarOx_;		//!< un-normalized mixture fraction on the oxidizer side

		std::vector<double> interpolated_Y_;
		std::vector<double> interpolated_omegap_;
		std::vector<double> interpolated_omegad_;

		std::vector<double> min_fields_Y_;	//!< minimum fields
		std::vector<double> max_fields_Y_;	//!< maximum fields
		std::vector<double> mean_fields_Y_;	//!< mean fields

		std::vector<double> min_fields_omegap_;			//!< minimum fields
		std::vector<double> max_fields_omegap_;			//!< maximum fields
		std::vector<double> mean_fields_omegap_;		//!< mean fields

		std::vector<double> min_fields_omegad_;		//!< minimum fields
		std::vector<double> max_fields_omegad_;		//!< maximum fields
		std::vector<double> mean_fields_omegad_;	//!< mean fields

		std::vector<bool> available_omega_;			//!< true if the production/consumption terms are available
		std::vector<unsigned int> list_available_omega_;	//!< list of available production/consumption terms 

		bool Y_def_PAH12_;		//!< true if PAH12 is available in the table
		bool Y_def_PAH34_;		//!< true if PAH34 is available in the table
		bool Y_def_PAHLP_;		//!< true if PAHLP is available in the table
		bool Y_def_SP_;			//!< true if SP is available in the table
		bool Y_def_AGG_;		//!< true if AGG is available in the table

		double Y_def_alpha_PAH12_;	//!< weight coefficient for PAH12
		double Y_def_alpha_PAH34_;	//!< weight coefficient for PAH34
		double Y_def_alpha_PAHLP_;	//!< weight coefficient for PAHLP
		double Y_def_alpha_SP_;		//!< weight coefficient for SP
		double Y_def_alpha_AGG_;	//!< weight coefficient for AGG

		std::vector<unsigned int> indices_pah12_;	//!< indices of species belonging to the PAH12 family
		std::vector<unsigned int> indices_pah34_;	//!< indices of species belonging to the PAH34 family
		std::vector<unsigned int> indices_pahlp_;	//!< indices of species belonging to the PAHLP family
		std::vector<unsigned int> indices_sp_;		//!< indices of species belonging to the SP family
		std::vector<unsigned int> indices_agg_;		//!< indices of species belonging to the AGG family
	};
}

#include "LookUpTable_ZY.hpp"

#endif /* OpenSMOKE_LookUpTable_ZY_H */

