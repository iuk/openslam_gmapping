#include <gmapping/gridfastslam/motionmodel.h>
#include <gmapping/utils/stat.h>
#include <iostream>

#define MotionModelConditioningLinearCovariance 0.01
#define MotionModelConditioningAngularCovariance 0.001

namespace GMapping {



OrientedPoint
MotionModel::drawFromMotion (const OrientedPoint& p, double linearMove, double angularMove) const{
	OrientedPoint n(p);
	double lm=linearMove  + fabs( linearMove ) * sampleGaussian( srr ) + fabs( angularMove ) * sampleGaussian( str );
	double am=angularMove + fabs( linearMove ) * sampleGaussian( srt ) + fabs( angularMove ) * sampleGaussian( stt );
	n.x+=lm*cos(n.theta+.5*am);
	n.y+=lm*sin(n.theta+.5*am);
	n.theta+=am;
	n.theta=atan2(sin(n.theta), cos(n.theta));
	return n;
}
// 输入的是：p：t-1 时刻的位姿（位置+角度），pnew：t 时刻的 odom 读数，pold t-1 时刻 odom 读数
// 返回的是：依据 odom 测量值 + 高斯噪音 + t-1 时刻位姿 得到的 t 时刻位姿
OrientedPoint MotionModel::drawFromMotion(
		const OrientedPoint& p,
		const OrientedPoint& pnew,
		const OrientedPoint& pold) const
{
	double sxy=0.3*srr;
	OrientedPoint delta=absoluteDifference(pnew, pold);
	OrientedPoint noisypoint(delta);	// 算出 odom 的增长值
	// fabs() 就是取绝对值
	// 下面操作 对 odom 中获取的偏移量 加入 高斯噪声
	// sampleGaussian 函数，按照 高斯分布来采样  在 stat.c 中，输入参数是方差，第二个参数是均值，默认 0
	noisypoint.x+=sampleGaussian(srr*fabs(delta.x)+str*fabs(delta.theta)+sxy*fabs(delta.y));
	noisypoint.y+=sampleGaussian(srr*fabs(delta.y)+str*fabs(delta.theta)+sxy*fabs(delta.x));
	noisypoint.theta+=sampleGaussian(stt*fabs(delta.theta)+srt*sqrt(delta.x*delta.x+delta.y*delta.y));	// 设 odom 值符合 高斯分布，给这个值加点符合高斯分布的随机数噪音干扰。

	// 加完噪音还是一组数阿，无法表达出这是高斯分布呀。还不如不加，不加的话的值应该是概率最大的值吧？？
	noisypoint.theta=fmod(noisypoint.theta, 2*M_PI);

	if (noisypoint.theta>M_PI)
		noisypoint.theta-=2*M_PI;

	// 返回的是依据 odom 测量值 + 高斯噪音 + t-1 时刻位姿 得到的 t 时刻位姿
	return absoluteSum(p,noisypoint);
}


/*
OrientedPoint
MotionModel::drawFromMotion(const OrientedPoint& p, const OrientedPoint& pnew, const OrientedPoint& pold) const{

	//compute the three stps needed for perfectly matching the two poses if the noise is absent

	OrientedPoint delta=pnew-pold;
	double aoffset=atan2(delta.y, delta.x);
	double alpha1=aoffset-pold.theta;
	alpha1=atan2(sin(alpha1), cos(alpha1));
	double rho=sqrt(delta*delta);
	double alpha2=pnew.theta-aoffset;
	alpha2=atan2(sin(alpha2), cos(alpha2));

	OrientedPoint pret=drawFromMotion(p, 0, alpha1);
	pret=drawFromMotion(pret, rho, 0);
	pret=drawFromMotion(pret, 0, alpha2);
	return pret;
}
*/


Covariance3 MotionModel::gaussianApproximation(const OrientedPoint& pnew, const OrientedPoint& pold) const{
	OrientedPoint delta=absoluteDifference(pnew,pold);
	double linearMove=sqrt(delta.x*delta.x+delta.y*delta.y);
	double angularMove=fabs(delta.x);
	double s11=srr*srr*linearMove*linearMove;
	double s22=stt*stt*angularMove*angularMove;
	double s12=str*angularMove*srt*linearMove;
	Covariance3 cov;
	double s=sin(pold.theta),c=cos(pold.theta);
	cov.xx=c*c*s11+MotionModelConditioningLinearCovariance;
	cov.yy=s*s*s11+MotionModelConditioningLinearCovariance;
	cov.tt=s22+MotionModelConditioningAngularCovariance;
	cov.xy=s*c*s11;
	cov.xt=c*s12;
	cov.yt=s*s12;
	return cov;
}

};
