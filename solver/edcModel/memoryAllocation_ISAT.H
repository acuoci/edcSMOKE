
Eigen::VectorXd phiISAT_HOM;
Eigen::VectorXd RphiISAT_HOM;
Eigen::MatrixXd mapGrad_HOM;

if(isatCheck == true)
{
	phiISAT_HOM.resize(thermodynamicsMapXML->NumberOfSpecies()+2);
	RphiISAT_HOM.resize(thermodynamicsMapXML->NumberOfSpecies()+2);
	mapGrad_HOM.resize(thermodynamicsMapXML->NumberOfSpecies()+2,thermodynamicsMapXML->NumberOfSpecies()+2);
}
