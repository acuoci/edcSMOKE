if (homogeneousReactions == true)
{
	if (combustionModel == "ED" || combustionModel == "EDFR")
		#include "EDFR.H"
	else if (combustionModel == "EDC")
		#include "EDC.H"
}
