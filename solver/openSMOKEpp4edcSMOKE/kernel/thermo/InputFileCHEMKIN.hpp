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

namespace OpenSMOKE
{
	bool ParseReactionClassLine(const std::string line_to_parse, const std::string keyword, std::vector<std::string>& names)
	{
		std::string line = line_to_parse;
		boost::replace_all(line, keyword, "");

		// Remove empty spaces
		boost::erase_all(line, " ");
		boost::erase_all(line, "\t");

		// Remove additional comments
		const size_t pos = line.find_last_of("]");
		if (pos != std::string::npos)
			line.erase(pos+1, line.length()-pos);

		// Remove special characters for comments
		boost::erase_all(line, "!");

		// Look for the existence of mainclass and subclass definitions
		if (line.find("][") == std::string::npos)
			return false;

		// Replace "][" with empty spaces
		boost::replace_all(line, "][", " ");

		// Remove square brackets
		boost::erase_all(line, "[");
		boost::erase_all(line, "]");

		// Separate mainclass and subclass
		boost::split(names, line, boost::is_any_of(" "));

		// Check if multiple definitions exist
		if (names.size() != 2)
			return false;

		return true;
	}

	InputFileCHEMKIN::InputFileCHEMKIN() 
	{
	}

	InputFileCHEMKIN::InputFileCHEMKIN(const std::string file_name) 
	{
            #if defined __linux || defined __APPLE__
                        
				//Check the file format through the awk command
				const std::string check_string = "if awk  '/\\r$/{exit 0;} 1{exit 1;}' " + file_name + 
					"; then exit 0; else  exit 1; fi; ";
				const int ff_result = system(check_string.c_str()) / 256;

				if (ff_result == 0)
				{
							// Check if the dos2unix application exists
							{
								const std::string exec_test = "which dos2unix >/dev/null";
								const int result_test = system(exec_test.c_str()) / 256;
					
								if (result_test == 0)
								{
									const std::string exec_dos2unix = "dos2unix " + file_name + " 2>/dev/null";
									const int result_dos2unix = system(exec_dos2unix.c_str());
								}
								else
								{
									const std::string exec_sed = "perl -pi -e 's/\r\n|\n|\r/\n/g' " + file_name;
									const int result_sed = system(exec_sed.c_str());
								}
							}
				}
                      
            #endif
			
            file_name_ = new boost::filesystem::path(file_name);
        
			std::ifstream myfile(file_name.c_str(), std::ios::in);
			CheckIfFileIsOpen(myfile, file_name);

			bool TakeFirstSpecies = false;
			int count=1;
			std::string line;
			while ( myfile.good() )
			{
					std::getline(myfile,line);

					size_t found_comment = line.find("!#");
					std::string comment = "";
					bool found_strong_comment = false;
					if (found_comment != line.npos)
					{
						comment = line.substr(found_comment);
						found_strong_comment = true;
					}

					size_t found_soot_class = line.find("[SOOTCLASS]");
					if (found_soot_class != line.npos)
					{
						std::string line_to_parse = line.substr(found_soot_class);
						boost::replace_all(line_to_parse, "[SOOTCLASS]", "");
						boost::erase_all(line_to_parse, " ");
						boost::erase_all(line_to_parse, "!");
						boost::erase_all(line_to_parse, "[");
						boost::erase_all(line_to_parse, "]");
						soot_class_names_.push_back(line_to_parse);
						soot_class_lines_abs_.push_back(count);
					}

					size_t found_reaction_class = line.find("[REACTIONCLASS]");
					if (found_reaction_class != line.npos)
					{
						std::string line_to_parse = line.substr(found_reaction_class);

						std::vector<std::string> classnames;
						const bool flag = ParseReactionClassLine(line_to_parse, "[REACTIONCLASS]", classnames);
						if (flag == false)
						{
							std::cout << "Line: " << count << line_to_parse << std::endl;
							OpenSMOKE::FatalErrorMessage("Syntax error in definition of [REACTIONCLASS]");
						}

						reaction_mainclass_names_.push_back(classnames[0]);
						reaction_mainclass_lines_abs_.push_back(count);
						reaction_subclass_names_.push_back(classnames[1]);
						reaction_subclass_lines_abs_.push_back(count);
					}

					size_t found_endreaction_class = line.find("[ENDREACTIONCLASS]");
					if (found_endreaction_class != line.npos)
					{
						std::string line_to_parse = line.substr(found_endreaction_class);

						std::vector<std::string> endclassnames;
						const bool flag = ParseReactionClassLine(line_to_parse, "[ENDREACTIONCLASS]", endclassnames);
						if (flag == false)
						{
							std::cout << "Line: " << count << line_to_parse << std::endl;
							OpenSMOKE::FatalErrorMessage("Syntax error in definition of [ENDREACTIONCLASS]");
						}

						endreaction_mainclass_names_.push_back(endclassnames[0]);
						endreaction_mainclass_lines_abs_.push_back(count);
						endreaction_subclass_names_.push_back(endclassnames[1]);
						endreaction_subclass_lines_abs_.push_back(count);
					}
					
					// AN to get the first species of each [SPECIESCLASS]
					if (TakeFirstSpecies == true) {
						boost::replace_all(line, "\t", " ");
						size_t found=line.find_first_of("!");
						if (found!=line.npos)
						line.erase(found);
						auto pos = line.find_first_not_of (' ');
						if (pos != line.npos) {				
							std::vector<std::string> stringvector;
							boost::split(stringvector,line.substr(pos),boost::is_any_of(" "));
							firstspecies_names_.push_back(stringvector.front());
							TakeFirstSpecies = false;
						}
					}

					size_t found_species_class = line.find("[SPECIESCLASS]");
					if (found_species_class != line.npos)
					{
						std::string line_to_parse = line.substr(found_species_class);
						boost::replace_all(line_to_parse, "[SPECIESCLASS]", "");
						boost::erase_all(line_to_parse, " ");
						boost::erase_all(line_to_parse, "!#");
						boost::erase_all(line_to_parse, "[");
						boost::erase_all(line_to_parse, "]");
						species_class_names_.push_back(line_to_parse);
						species_class_lines_abs_.push_back(count);
						//AN to get the first species of each [SPECIESCLASS]
						TakeFirstSpecies = true;
					}

					size_t found=line.find_first_of("!");
					if (found!=line.npos)
						line.erase(found);

					// replacing tabs with spaces
					boost::replace_all(line, "\t", " ");
					if (line.find_first_not_of (' ') == line.npos)   // has only spaces?
					{
						indices_of_blank_lines_.push_back(count);

						if (found_strong_comment == true)
						{
							indices_of_strong_comment_lines_.push_back(count);
							strong_comment_lines_.push_back(comment);
						}
					}
					else
					{	
						good_lines_.push_back(line);
						indices_of_good_lines_.push_back(count);
						strong_comments_.push_back(comment);
					}
					count++;
			}
			myfile.close();

        
			number_of_blank_lines_ = boost::lexical_cast<int>(indices_of_blank_lines_.size());
			number_of_good_lines_ = boost::lexical_cast<int>(indices_of_good_lines_.size());
			number_of_strong_comment_lines_ = boost::lexical_cast<int>(indices_of_strong_comment_lines_.size());
			number_of_lines_ = number_of_blank_lines_ + number_of_good_lines_;
	}

	void InputFileCHEMKIN::ConvertGoodLinesIntoBlankLines(const std::vector<unsigned int> lines_to_remove)
	{
		std::vector<unsigned int> indices = lines_to_remove;
		std::sort(indices.begin(), indices.end());
		std::reverse(indices.begin(), indices.end());

		number_of_blank_lines_ += static_cast<int>(indices.size());
		number_of_good_lines_ -= static_cast<int>(indices.size());

		for (unsigned int j = 0; j < indices.size(); j++)
		{
			indices_of_blank_lines_.push_back(indices_of_good_lines_[indices[j]]);
			indices_of_good_lines_.erase(indices_of_good_lines_.begin() + (indices[j]-1));
		}

		for (unsigned int j = 0; j<indices.size(); j++)
			good_lines_.erase(good_lines_.begin() + (indices[j]-1));	
	}

	InputFileCHEMKIN::InputFileCHEMKIN(const InputFileCHEMKIN& orig) {
	}

	InputFileCHEMKIN::~InputFileCHEMKIN() {
	}

	void InputFileCHEMKIN::Status(std::ostream &fOut) const
	{ 
		fOut << "Name:        " << file_name_->filename() << std::endl;
		fOut << "Path:        " << file_name_->parent_path() << std::endl;
		fOut << "Size:        " << boost::filesystem::file_size(*file_name_)/1000. << " kB" << std::endl;
		fOut << "Lines:       " << number_of_lines_ << std::endl;
		fOut << "Blank lines: " << number_of_blank_lines_ << std::endl;
	}

}
