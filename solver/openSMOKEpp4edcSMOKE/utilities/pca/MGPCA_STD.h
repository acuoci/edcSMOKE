/*----------------------------------------------------------------------*\
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
|   Copyright(C) 2014, 2013, 2012  Alberto Cuoci                          |
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

#ifndef OpenSMOKE_MGPCA_STD_H
#define OpenSMOKE_MGPCA_STD_H

#include "math/OpenSMOKEClass.hpp"
#include "math/OpenSMOKEVector.h"

namespace OpenSMOKE
{
	//!  A class to perform MGPCA Standard (without clustering)
	/*!
			A class to perform MGPCA Standard (without clustering)
	*/

	template <typename ThermoMap>
	class MGPCA_STD
	{
	
	public:

		enum MGPCA_Type { MGPCA_Standard, MGPCA_Standard_LocalScaling };

		MGPCA_STD(const ThermoMap& thermo, const boost::filesystem::path& folder_name, const MGPCA_Type type);
	
		void SetNormalization(const bool normalization);
		
		void Reconstruct(Eigen::VectorXd& x);
		void Reconstruct(const double z, Eigen::VectorXd& x);

		MGPCA_Type type() const { return type_; };

		unsigned int number_of_retained_species() const { return number_of_retained_species_; }
		unsigned int number_of_unretained_species() const { return number_of_unretained_species_; }
		unsigned int number_of_removed_species() const { return number_of_removed_species_; }
		unsigned int number_of_species() const { return number_of_species_; }

		const std::vector<unsigned int>& indices_of_retained_species() const { return indices_of_retained_species_; }
		const std::vector<unsigned int>& indices_of_unretained_species() const { return indices_of_unretained_species_; }
		const std::vector<unsigned int>& indices_of_removed_species() const { return indices_of_removed_species_; }

	protected:

		unsigned int number_of_retained_species_;		//!< number of retained species (PV)
		unsigned int number_of_unretained_species_;		//!< number of unretained species
		unsigned int number_of_removed_species_;		//!< number of removed species
		unsigned int number_of_tracked_species_;		//!< total number of tracked species
		unsigned int number_of_species_;				//!< total number of species

		std::vector<unsigned int> indices_of_retained_species_;		//!< indices of retained species (0-based)
		std::vector<unsigned int> indices_of_unretained_species_;	//!< indices of unretained species (0-based)
		std::vector<unsigned int> indices_of_removed_species_;		//!< indices of removed species (0-based)
		std::vector<unsigned int> indices_of_tracked_species_;		//!< indices of tracked species (0-based)

		std::vector< Eigen::VectorXd > average_retained_;			//!< vector of averages of retained species
		std::vector< Eigen::VectorXd > average_unretained_;			//!< vector of averages of retained species

		Eigen::VectorXd gamma_retained_;			//!< vector of scaling factors of retained species
		Eigen::VectorXd gamma_unretained_;			//!< vector of scaling factors of retained species
		
		Eigen::MatrixXd B_;						//!< projection matrix
		Eigen::VectorXd x_retained_;			//!< vector of retained species
		Eigen::VectorXd x_unretained_;			//!< vector of unretained species

		bool normalization_;
		MGPCA_Type type_;
		unsigned int number_of_mixture_fraction_points_;

		Eigen::VectorXd current_average_retained_;
		Eigen::VectorXd current_average_unretained_;
	};
}

#include "MGPCA_STD.hpp"

#endif // OpenSMOKE_MGPCA_STD_H
