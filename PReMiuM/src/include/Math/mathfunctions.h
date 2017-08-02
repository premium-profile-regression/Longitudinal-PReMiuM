/// \file mathfunctions.h
/// \author David Hastie
/// \date 19 Mar 2012
/// \brief Header file to define distributions

/// \note (C) Copyright David Hastie and Silvia Liverani, 2012.

/// PReMiuM++ is free software; you can redistribute it and/or modify it under the
/// terms of the GNU Lesser General Public License as published by the Free Software
/// Foundation; either version 3 of the License, or (at your option) any later
/// version.

/// PReMiuM++ is distributed in the hope that it will be useful, but WITHOUT ANY
/// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
/// PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

/// You should have received a copy of the GNU Lesser General Public License
/// along with PReMiuM++ in the documentation directory. If not, see
/// <http://www.gnu.org/licenses/>.

/// The external linear algebra library Eigen, parts of which are included  in the
/// lib directory is released under the LGPL3+ licence. See comments in file headers
/// for details.

/// The Boost C++ header library, parts of which are included in the  lib directory
/// is released under the Boost Software Licence, Version 1.0, a copy  of which is
/// included in the documentation directory.

/// Version 3.1.3 edited by Rob Johnson, March 2017


#ifndef MATHFUNCTIONS_H_
#define MATHFUNCTIONS_H_

#include<cmath>

#include<boost/math/special_functions/gamma.hpp>

using namespace boost::math::constants;

using boost::math::lgamma;

double logMultivarGammaFn(const double& x,const unsigned int& p){

	double out;
	out = 0.25*(double)(p*(p-1))*log(pi<double>());
	for(unsigned int i=1;i<=p;i++){
		out += lgamma(x+(1.0-(double)i)/2.0);
	}
	return out;
}

double logit(const double& lambda){
	return 1.0/(1.0+exp(-lambda));
}

//RJ matrix (fast) inversion function
 void invert(MatrixXd& invSigma,MatrixXd& Sigma,const unsigned int dimBlock, double noise){
	if(dimBlock>0){
		int dimSigma = Sigma.rows();
		int nBlocks = dimSigma/dimBlock;
		MatrixXd C;
		C = Sigma.block(0, 0, dimBlock, dimBlock);
		for(int i=0;i<dimBlock;i++)
			C(i,i) = C(i,i) - noise;
		double invNoise = 1.0/noise;
		MatrixXd E;
        E = nBlocks*C;
        for(int i=0;i<dimBlock;i++)
            E(i,i) = E(i,i) + noise;
		E = - invNoise * C * E.inverse();
        invSigma = E.replicate(nBlocks,nBlocks);
        for(int i=0;i<dimSigma;i++)
            invSigma(i,i) = invSigma(i,i) + invNoise;
	}else{
		invSigma = Sigma.inverse();
	}
}
//RJ GP covariance function
void GP_cov(MatrixXd& Mat,std::vector<double> L,std::vector<double> times,const unsigned int dimBlock){
	double a;
  	int i,j;
  	int nTimes = times.size();
  	double eL0 = exp(L[0]);
  	double eL1 = exp(L[1])*2.0;
  	double eL2 = exp(L[2]);
  	if(dimBlock<0){
  		int nBlocks = nTimes/dimBlock;
		MatrixXd unit;
		unit.resize(dimBlock,dimBlock);
        unit.fill(0.0);
  		for(i=1;i<dimBlock;i++){
    		for(j=0;j<i;j++){
        		a=-(times[i]-times[j])*(times[i]-times[j])/eL1;
        		unit(i,j)=eL0*std::exp(a);
    		}
  		}
        unit = unit + unit.transpose() + eL0*MatrixXd::Identity(dimBlock, dimBlock);
  		Mat = unit.replicate(nBlocks,nBlocks);
        for(i=0;i<nTimes;i++) Mat(i,i) = Mat(i,i) + eL2;
    }else{
  	
		for(i=1;i<nTimes;i++){
    		for(j=0;j<i;j++){
        		a=-(times[i]-times[j])*(times[i]-times[j])/eL1;
        		Mat(i,j)=eL0*std::exp(a);
    		}
  		}	
        Mat = Mat + Mat.transpose();
        for(int i=0;i<nTimes;i++)
            Mat(i,i) = Mat(i,i) + eL0+eL2;
  	}
}
//RJ GP covariance star function
double GP_cov_star(MatrixXd& Mat,std::vector<double> L,std::vector<double> times,std::vector<double> timestar){
	
	double a;
  	int i,j;
  	int nTimes = times.size();
  	int nStar = timestar.size();
  	for(i=0;i<nTimes;i++){
    	for(j=0;j<nStar;j++){
        	a=-1.0/2.0/exp(L[1])*(times[i]-timestar[j])*(times[i]-timestar[j]);
        	Mat(i,j)=exp(L[0])*std::exp(a);
    	}
  	}	
  	return 0;
}

#endif /*MATHFUNCTIONS_H_*/
