#ifndef ODESystemVirtualClassWithOpenSMOKEVectors_H
#define	ODESystemVirtualClassWithOpenSMOKEVectors_H

namespace OpenSMOKE
{
	
	class ODESystemVirtualClassWithOpenSMOKEVectors
	{
	public:

		unsigned int NumberOfEquations() { return ne_; }
		virtual void GetEquations(const OpenSMOKE::OpenSMOKEVectorDouble& y, const double t, OpenSMOKE::OpenSMOKEVectorDouble& dy) = 0;
		virtual void PrintResults(const OpenSMOKE::OpenSMOKEVectorDouble& y, const double t) { };

	protected:

		unsigned int ne_;

		void MemoryAllocation()
		{
			ChangeDimensions(ne_, &y_, true);
			ChangeDimensions(ne_, &dy_, false);
		}

		void Equations(const Eigen::VectorXd &Y, const double t, Eigen::VectorXd &DY)
		{
			y_.CopyFrom(Y.data());
			GetEquations(y_, t, dy_);
			dy_.CopyTo(DY.data());
		}

		void Jacobian(const Eigen::VectorXd &Y, const double t, Eigen::MatrixXd &J) { };

		void Print(const double t, const Eigen::VectorXd &Y)
		{
			y_.CopyFrom(Y.data());
			PrintResults(y_, t);
		}

	private:

		OpenSMOKE::OpenSMOKEVectorDouble  y_;
		OpenSMOKE::OpenSMOKEVectorDouble dy_;
	};
}
	
#endif	// ODESystemVirtualClassWithOpenSMOKEVectors
