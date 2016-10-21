
Eigen::VectorXd phiISAT_HOM;
Eigen::VectorXd RphiISAT_HOM;
Eigen::MatrixXd mapGrad_HOM;

if(isatCheck == true)
{
	phiISAT_HOM.resize(ne_total);
	RphiISAT_HOM.resize(ne_total);
	mapGrad_HOM.resize(ne_total,ne_total);
}
