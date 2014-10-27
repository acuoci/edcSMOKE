
//- HOM
Eigen::VectorXd phiISAT_HOM;
Eigen::VectorXd RphiISAT_HOM;
Eigen::MatrixXd mapGrad_HOM;

//- Stats file
int counterISATHOM = 0;
FILE *fleHOM, *fleTme, *fleCht, *fleChm;

if(isatCheck == true)
{
	phiISAT_HOM.resize(thermodynamicsMapXML->NumberOfSpecies()+2);
	RphiISAT_HOM.resize(thermodynamicsMapXML->NumberOfSpecies()+2);
	mapGrad_HOM.resize(thermodynamicsMapXML->NumberOfSpecies()+2,thermodynamicsMapXML->NumberOfSpecies()+2);

	if(printStats == true) 
	{
		ostringstream st;
		st << getpid();
		std::string strtme = "tmeExec_";
		strtme.append(st.str());
		fleTme = fopen(strtme.c_str(),"w+");

		{
			std::string strchm = "tmeChem_HOM_";
			strchm.append(st.str());
			fleChm = fopen(strchm.c_str(),"w+");
		}
	
		if(isatCheck == true) 
		{
			{
				std::string strhom = "statsHOM_";
				strhom.append(st.str());
				fleHOM = fopen(strhom.c_str(),"w+");
			}
		}
	}
}
