  #include <cstring>
#include <limits>
#include <list>
#include <iostream>

#include <gmapping/scanmatcher/scanmatcher.h>
#include "gridlinetraversal.h"
//#define GENERATE_MAPS

namespace GMapping {

using namespace std;

const double ScanMatcher::nullLikelihood=-.5;

ScanMatcher::ScanMatcher(): m_laserPose(0,0,0){
	//m_laserAngles=0;
	m_laserBeams=0;
	m_optRecursiveIterations=3;
	m_activeAreaComputed=false;

	// This  are the dafault settings for a grid map of 5 cm
	m_llsamplerange=0.01;
	m_llsamplestep=0.01;
	m_lasamplerange=0.005;
	m_lasamplestep=0.005;
	m_enlargeStep=10.;
	m_fullnessThreshold=0.1;
    m_angularOdometryReliability=0.;    //
    m_linearOdometryReliability=0.;     //
    m_freeCellRatio=sqrt(2.);   // 对角线的意思把
	m_initialBeamsSkip=0;

/*
	// This  are the dafault settings for a grid map of 10 cm
	m_llsamplerange=0.1;
	m_llsamplestep=0.1;
	m_lasamplerange=0.02;
	m_lasamplestep=0.01;
*/
	// This  are the dafault settings for a grid map of 20/25 cm
/*
	m_llsamplerange=0.2;
	m_llsamplestep=0.1;
	m_lasamplerange=0.02;
	m_lasamplestep=0.01;
	m_generateMap=false;
*/

   m_linePoints = new IntPoint[20000];
}

ScanMatcher::~ScanMatcher(){
	delete [] m_linePoints;
}

void ScanMatcher::invalidateActiveArea(){
	m_activeAreaComputed=false;
}

/*
void ScanMatcher::computeActiveArea(ScanMatcherMap& map, const OrientedPoint& p, const double* readings){
	if (m_activeAreaComputed)
		return;
	HierarchicalArray2D<PointAccumulator>::PointSet activeArea;
	OrientedPoint lp=p;
	lp.x+=cos(p.theta)*m_laserPose.x-sin(p.theta)*m_laserPose.y;
	lp.y+=sin(p.theta)*m_laserPose.x+cos(p.theta)*m_laserPose.y;
	lp.theta+=m_laserPose.theta;
	IntPoint p0=map.world2map(lp);
	const double * angle=m_laserAngles;
	for (const double* r=readings; r<readings+m_laserBeams; r++, angle++)
		if (m_generateMap){
			double d=*r;
			if (d>m_laserMaxRange)
				continue;
			if (d>m_usableRange)
				d=m_usableRange;

			Point phit=lp+Point(d*cos(lp.theta+*angle),d*sin(lp.theta+*angle));
			IntPoint p1=map.world2map(phit);

			d+=map.getDelta();
			//Point phit2=lp+Point(d*cos(lp.theta+*angle),d*sin(lp.theta+*angle));
			//IntPoint p2=map.world2map(phit2);
			IntPoint linePoints[20000] ;
			GridLineTraversalLine line;
			line.points=linePoints;
			//GridLineTraversal::gridLine(p0, p2, &line);
			GridLineTraversal::gridLine(p0, p1, &line);
			for (int i=0; i<line.num_points-1; i++){
				activeArea.insert(map.storage().patchIndexes(linePoints[i]));
			}
			if (d<=m_usableRange){
				activeArea.insert(map.storage().patchIndexes(p1));
				//activeArea.insert(map.storage().patchIndexes(p2));
			}
		} else {
			if (*r>m_laserMaxRange||*r>m_usableRange) continue;
			Point phit=lp;
			phit.x+=*r*cos(lp.theta+*angle);
			phit.y+=*r*sin(lp.theta+*angle);
			IntPoint p1=map.world2map(phit);
			assert(p1.x>=0 && p1.y>=0);
			IntPoint cp=map.storage().patchIndexes(p1);
			assert(cp.x>=0 && cp.y>=0);
			activeArea.insert(cp);

		}
	//this allocates the unallocated cells in the active area of the map
	//cout << "activeArea::size() " << activeArea.size() << endl;
	map.storage().setActiveArea(activeArea, true);
	m_activeAreaComputed=true;
}
*/

// 计算有效区域
// 设置 粒子的地图，只是当前粒子 当前扫描得到的地图，不是全部历史的地图？？
// ScanMatcherMap ： 粒子中内存的地图
// OrientedPoint：粒子当前的位置（已被更新为 t 时刻）
// readings：t 时刻的扫描数据
void ScanMatcher::computeActiveArea(ScanMatcherMap& map,
                                    const OrientedPoint& p, // 粒子的位置，在世界中的位置
                                    const double* readings)
{
	if (m_activeAreaComputed)	// 标志位，如果已经被计算过了
		return;

	OrientedPoint lp=p;	// laser 在 world 中的位置
	//m_laserPose 就是 laser 与 base_link 的静态相对位置吧。
	lp.x+=cos(p.theta)*m_laserPose.x-sin(p.theta)*m_laserPose.y;
	lp.y+=sin(p.theta)*m_laserPose.x+cos(p.theta)*m_laserPose.y;
	lp.theta+=m_laserPose.theta;
    IntPoint p0=map.world2map(lp);  // p0，laser中心在 map 中的位置

    Point min(map.map2world(0,0));  // world 的 最小坐标角的坐标
    Point max(map.map2world(map.getMapSizeX()-1,map.getMapSizeY()-1));  // world 中的最大坐标

    if (lp.x<min.x) min.x=lp.x; // 根据 laser 位姿，看看是否需要修改 world 的最小最大点
	if (lp.y<min.y) min.y=lp.y;
	if (lp.x>max.x) max.x=lp.x;
	if (lp.y>max.y) max.y=lp.y;

	/*determine the size of the area*/
	const double * angle=m_laserAngles+m_initialBeamsSkip;

    // 遍历所有激光点
    // 目的是：修改 world 的最大最小点？
    for (const double* r=readings+m_initialBeamsSkip; r<readings+m_laserBeams; r++, angle++)
    {
        if (*r>m_laserMaxRange||*r==0.0||isnan(*r)) continue;   // 排除异常激光点

        double d=*r>m_usableRange?m_usableRange:*r;
        //由上可见 大于 MaxRange 被排除，大于 usableRange 被 =usableRange

        Point phit=lp;  // 计算激光点的坐标，在世界中的
		phit.x+=d*cos(lp.theta+*angle);
		phit.y+=d*sin(lp.theta+*angle);
        if (phit.x<min.x) min.x=phit.x; // 根据激光点位置，扩张 world 的最大最小点
		if (phit.y<min.y) min.y=phit.y;
		if (phit.x>max.x) max.x=phit.x;
		if (phit.y>max.y) max.y=phit.y;
	}
	//min=min-Point(map.getDelta(),map.getDelta());
	//max=max+Point(map.getDelta(),map.getDelta());

    if ( !map.isInside(min)	|| !map.isInside(max))  // 如果 地图 有扩张
	{
		Point lmin(map.map2world(0,0));
		Point lmax(map.map2world(map.getMapSizeX()-1,map.getMapSizeY()-1));
		//cerr << "CURRENT MAP " << lmin.x << " " << lmin.y << " " << lmax.x << " " << lmax.y << endl;
		//cerr << "BOUNDARY OVERRIDE " << min.x << " " << min.y << " " << max.x << " " << max.y << endl;
		min.x=( min.x >= lmin.x )? lmin.x: min.x-m_enlargeStep;
		max.x=( max.x <= lmax.x )? lmax.x: max.x+m_enlargeStep;
		min.y=( min.y >= lmin.y )? lmin.y: min.y-m_enlargeStep;
		max.y=( max.y <= lmax.y )? lmax.y: max.y+m_enlargeStep;
        map.resize(min.x, min.y, max.x, max.y);     // 地图 重定义大小
		//cerr << "RESIZE " << min.x << " " << min.y << " " << max.x << " " << max.y << endl;
	}

    HierarchicalArray2D<PointAccumulator>::PointSet activeArea; // 这是地图所用的储存结构，同款
	/*allocate the active area*/
	angle=m_laserAngles+m_initialBeamsSkip;

    // 又一次 遍历所有 激光点 数据
	for (const double* r=readings+m_initialBeamsSkip; r<readings+m_laserBeams; r++, angle++)
	{
		if (m_generateMap)	// bool 如果 生成地图
		{
			double d=*r;	// 激光点数据就是距离数据，存到 d 中

			if (d>m_laserMaxRange||d==0.0||isnan(d))	// 如果这个 d 异常，则 跳过
				continue;

			if (d>m_usableRange)	// 如果 d 超过 有效距离最大值，则令其=有效距离最大值
				d=m_usableRange;

            // 激光点的坐标  // 把距离转化为 激光点的坐标（world中）
            Point phit=lp+Point(d*cos(lp.theta+*angle),d*sin(lp.theta+*angle));

			IntPoint p0=map.world2map(lp);	// laser 在 map 中的坐标
            IntPoint p1=map.world2map(phit);// 激光点 在 map 中的坐标

			//IntPoint linePoints[20000] ;
			GridLineTraversalLine line;
            line.points=m_linePoints;   // IntPoint 的 vector
			GridLineTraversal::gridLine(p0, p1, &line);	// 计算两点连线中的点？？存放在 line 中

			for (int i=0; i<line.num_points-1; i++)
			{
				assert(map.isInside(m_linePoints[i]));
                // 把两点连线间的点插入到 新建的 activeArea 地图 中
                // patchIndexes 输入 point(float型)，输出 IntPoint
				activeArea.insert(map.storage().patchIndexes(m_linePoints[i]));
				assert(m_linePoints[i].x>=0 && m_linePoints[i].y>=0);
			}

			if (d<m_usableRange)
			{
				IntPoint cp=map.storage().patchIndexes(p1);
				assert(cp.x>=0 && cp.y>=0);
				activeArea.insert(cp);
			}
		} else	// 如果 不生成 地图
		{
			if (*r>m_laserMaxRange||*r>m_usableRange||*r==0.0||isnan(*r)) continue;
			Point phit=lp;
			phit.x+=*r*cos(lp.theta+*angle);
			phit.y+=*r*sin(lp.theta+*angle);
			IntPoint p1=map.world2map(phit);
			assert(p1.x>=0 && p1.y>=0);
			IntPoint cp=map.storage().patchIndexes(p1);
			assert(cp.x>=0 && cp.y>=0);
			activeArea.insert(cp);
		}
	}

	//this allocates the unallocated cells in the active area of the map
	//cout << "activeArea::size() " << activeArea.size() << endl;
/*
	cerr << "ActiveArea=";
	for (HierarchicalArray2D<PointAccumulator>::PointSet::const_iterator it=activeArea.begin(); it!= activeArea.end(); it++){
		cerr << "(" << it->x <<"," << it->y << ") ";
	}
	cerr << endl;
*/

// setActiveArea 位于 harry2d.h 中
	map.storage().setActiveArea(activeArea, true);	// 设置 粒子的地图，只是当前粒子 当前扫描得到的地图，不是全部历史的地图？
	m_activeAreaComputed=true;
}

// 计算 地图每个点的 占用概率
// ScanMatcherMap ： 粒子中内存的地图
// OrientedPoint：粒子当前的位置（已被更新为 t 时刻）
// readings：t 时刻的扫描数据D
double ScanMatcher::registerScan(ScanMatcherMap& map, const OrientedPoint& p, const double* readings)
{
	if (!m_activeAreaComputed)
		computeActiveArea(map, p, readings);

	//this operation replicates the cells that will be changed in the registration operation
	map.storage().allocActiveArea();

	OrientedPoint lp=p;
	lp.x+=cos(p.theta)*m_laserPose.x-sin(p.theta)*m_laserPose.y;
	lp.y+=sin(p.theta)*m_laserPose.x+cos(p.theta)*m_laserPose.y;
	lp.theta+=m_laserPose.theta;
	IntPoint p0=map.world2map(lp);


	const double * angle=m_laserAngles+m_initialBeamsSkip;
	double esum=0;

	// 遍历所有激光数据点
	for (const double* r=readings+m_initialBeamsSkip; r<readings+m_laserBeams; r++, angle++)
	{
		if (m_generateMap)
		{
			double d=*r;
			if (d>m_laserMaxRange||d==0.0||isnan(d))
				continue;
			if (d>m_usableRange)
				d=m_usableRange;
			Point phit=lp+Point(d*cos(lp.theta+*angle),d*sin(lp.theta+*angle));
			IntPoint p1=map.world2map(phit);
			//IntPoint linePoints[20000] ;
			GridLineTraversalLine line;
			line.points=m_linePoints;
			GridLineTraversal::gridLine(p0, p1, &line);
			for (int i=0; i<line.num_points-1; i++)
			{
				PointAccumulator& cell=map.cell(line.points[i]);
				double e=-cell.entropy();
				cell.update(false, Point(0,0));
				e+=cell.entropy();
				esum+=e;
			}
			if (d<m_usableRange)
			{
				double e=-map.cell(p1).entropy();
				map.cell(p1).update(true, phit);
				e+=map.cell(p1).entropy();
				esum+=e;
			}
		} else
		{
			if (*r>m_laserMaxRange||*r>m_usableRange||*r==0.0||isnan(*r)) continue;
			Point phit=lp;
			phit.x+=*r*cos(lp.theta+*angle);
			phit.y+=*r*sin(lp.theta+*angle);
			IntPoint p1=map.world2map(phit);
			assert(p1.x>=0 && p1.y>=0);
			map.cell(p1).update(true,phit);
		}
	}
	//cout  << "informationGain=" << -esum << endl;
	return esum;
}

/*
void ScanMatcher::registerScan(ScanMatcherMap& map, const OrientedPoint& p, const double* readings){
	if (!m_activeAreaComputed)
		computeActiveArea(map, p, readings);

	//this operation replicates the cells that will be changed in the registration operation
	map.storage().allocActiveArea();

	OrientedPoint lp=p;
	lp.x+=cos(p.theta)*m_laserPose.x-sin(p.theta)*m_laserPose.y;
	lp.y+=sin(p.theta)*m_laserPose.x+cos(p.theta)*m_laserPose.y;
	lp.theta+=m_laserPose.theta;
	IntPoint p0=map.world2map(lp);
	const double * angle=m_laserAngles;
	for (const double* r=readings; r<readings+m_laserBeams; r++, angle++)
		if (m_generateMap){
			double d=*r;
			if (d>m_laserMaxRange)
				continue;
			if (d>m_usableRange)
				d=m_usableRange;
			Point phit=lp+Point(d*cos(lp.theta+*angle),d*sin(lp.theta+*angle));
			IntPoint p1=map.world2map(phit);

			IntPoint linePoints[20000] ;
			GridLineTraversalLine line;
			line.points=linePoints;
			GridLineTraversal::gridLine(p0, p1, &line);
			for (int i=0; i<line.num_points-1; i++){
				IntPoint ci=map.storage().patchIndexes(line.points[i]);
				if (map.storage().getActiveArea().find(ci)==map.storage().getActiveArea().end())
					cerr << "BIG ERROR" <<endl;
				map.cell(line.points[i]).update(false, Point(0,0));
			}
			if (d<=m_usableRange){

				map.cell(p1).update(true,phit);
			}
		} else {
			if (*r>m_laserMaxRange||*r>m_usableRange) continue;
			Point phit=lp;
			phit.x+=*r*cos(lp.theta+*angle);
			phit.y+=*r*sin(lp.theta+*angle);
			map.cell(phit).update(true,phit);
		}
}

*/

double ScanMatcher::icpOptimize(OrientedPoint& pnew, const ScanMatcherMap& map, const OrientedPoint& init, const double* readings) const{
	double currentScore;
	double sc=score(map, init, readings);;
	OrientedPoint start=init;
	pnew=init;
	int iterations=0;
	do{
		currentScore=sc;
		sc=icpStep(pnew, map, start, readings);
		//cerr << "pstart=" << start.x << " " <<start.y << " " << start.theta << endl;
		//cerr << "pret=" << pnew.x << " " <<pnew.y << " " << pnew.theta << endl;
		start=pnew;
		iterations++;
	} while (sc>currentScore);
	cerr << "i="<< iterations << endl;
	return currentScore;
}

//输入：map 粒子携带地图，init t 时刻粒子位置（初始值），t 时刻 雷达数据
//输出：pnew 优化后的 t 时刻粒子位置，粒子得分 score
// 使用了一种简单迭代的方式，在初始 init 周围，寻找分手最高的 点
double ScanMatcher::optimize(OrientedPoint& pnew, const ScanMatcherMap& map,
                             const OrientedPoint& init, const double* readings) const
{
	double bestScore=-1;
	OrientedPoint currentPose=init;	// t 时刻 粒子的位置，作为初始值。这个位置是通过 odom 数据得来的。

    // 上来就计算一下初始位置的 score，score 的函数定义在 .h 文件中。
    // score 越小，则说明 scan 数据 与 map 越贴近
    // 但是 score 不会小于 0
    double currentScore=score(map, currentPose, readings);  // 重要函数

	double adelta=m_optAngularDelta, ldelta=m_optLinearDelta;
	unsigned int refinement=0;
	enum Move{Front, Back, Left, Right, TurnLeft, TurnRight, Done};
/*	cout << __PRETTY_FUNCTION__<<  " readings: ";
	for (int i=0; i<m_laserBeams; i++){
		cout << readings[i] << " ";
	}
	cout << endl;
*/	int c_iterations=0;

    do{ // 这个do 的：while (currentScore>bestScore || refinement<m_optRecursiveIterations);

		// 一开始肯定执行不了这句
		if (bestScore>=currentScore)	// 如果搜到的 currentscroe 比较小
		{
			refinement++;	// 搜小搜索范围，这是缩小次数的计数
			adelta*=.5;
			ldelta*=.5;
		}
//！！！！！ 是不是少了一个 else ？
		//一开始和以后都可以执行这句

		bestScore=currentScore;

//		cout <<"score="<< currentScore << " refinement=" << refinement;
//		cout <<  "pose=" << currentPose.x  << " " << currentPose.y << " " << currentPose.theta << endl;
		OrientedPoint bestLocalPose=currentPose;
		OrientedPoint localPose=currentPose;

		Move move=Front;
        do{ // 对应的：while(move!=Done);

			localPose=currentPose;

			// 把 粒子初始的 current pose 前后左右 平移，顺时 逆时 旋转，得到不同位置
			// 看看略微变动后，那个位姿的 得分 最 高

			switch(move)
			{
				case Front:
					localPose.x+=ldelta;
					move=Back;
					break;
				case Back:
					localPose.x-=ldelta;
					move=Left;
					break;
				case Left:
					localPose.y-=ldelta;
					move=Right;
					break;
				case Right:
					localPose.y+=ldelta;
					move=TurnLeft;
					break;
				case TurnLeft:
					localPose.theta+=adelta;
					move=TurnRight;
					break;
				case TurnRight:
					localPose.theta-=adelta;
					move=Done;
					break;
				default:;
            }
            // odo_gain 是什么..
			double odo_gain=1;
			if (m_angularOdometryReliability>0.)
			{
                double dth=init.theta-localPose.theta;  // 通过小步挪动，转了多少
                dth=atan2(sin(dth), cos(dth));  // 这是噶哈？
                dth*=dth;   // 自乘，所以一定>0
                // dth 越大，odo_gain 越小，且odo_gain总是（0，1）之间
                odo_gain*=exp(-m_angularOdometryReliability*dth);
			}
			if (m_linearOdometryReliability>0.)
			{
				double dx=init.x-localPose.x;
				double dy=init.y-localPose.y;
                double drho=dx*dx+dy*dy;    // 通过小步挪动，动了多少直线距离
                // drho 越大，odo_gain 越小，且odo_gain总是（0，1）之间
                odo_gain*=exp(-m_linearOdometryReliability*drho);
			}

			// 又做一次 score ?
            // odo_gain 的作用就是，如果挪动的距离远了，借给你强制降分，将分意味着相似度高呀？
            double localScore=odo_gain*score(map, localPose, readings);

            // currentscore 只取 localscore 的最大值
			if (localScore>currentScore)
			{
				currentScore=localScore;	// 最佳周围得分
				bestLocalPose=localPose;	// 最佳周围点
			}
			c_iterations++;
		} while(move!=Done);

		currentPose=bestLocalPose;
//		cout << "currentScore=" << currentScore<< endl;
		//here we look for the best move;
		// 构造函数中  m_optRecursiveIterations = 3
		// 迭代次数 要 >= 3 次
	}while (currentScore>bestScore || refinement<m_optRecursiveIterations);
    // 通过一次取周围样本点，得到一个分最高的，再在这个高分点周围取点迭代。
    // 终止条件是：1.周围点的相似度都小于中心点。并且 2.已经缩小三次搜索范围了
	//cout << __PRETTY_FUNCTION__ << "bestScore=" << bestScore<< endl;
	//cout << __PRETTY_FUNCTION__ << "iterations=" << c_iterations<< endl;

	pnew=currentPose;
	return bestScore;
}

struct ScoredMove{
	OrientedPoint pose;
	double score;
	double likelihood;
};

typedef std::list<ScoredMove> ScoredMoveList;

double ScanMatcher::optimize(OrientedPoint& _mean, ScanMatcher::CovarianceMatrix& _cov, const ScanMatcherMap& map, const OrientedPoint& init, const double* readings) const{
	ScoredMoveList moveList;
	double bestScore=-1;
	OrientedPoint currentPose=init;
	ScoredMove sm={currentPose,0,0};
	unsigned int matched=likelihoodAndScore(sm.score, sm.likelihood, map, currentPose, readings);
	double currentScore=sm.score;
	moveList.push_back(sm);
	double adelta=m_optAngularDelta, ldelta=m_optLinearDelta;
	unsigned int refinement=0;
	int count=0;
	enum Move{Front, Back, Left, Right, TurnLeft, TurnRight, Done};
	do{
		if (bestScore>=currentScore){
			refinement++;
			adelta*=.5;
			ldelta*=.5;
		}
		bestScore=currentScore;
//		cout <<"score="<< currentScore << " refinement=" << refinement;
//		cout <<  "pose=" << currentPose.x  << " " << currentPose.y << " " << currentPose.theta << endl;
		OrientedPoint bestLocalPose=currentPose;
		OrientedPoint localPose=currentPose;

		Move move=Front;
		do {
			localPose=currentPose;
			switch(move){
				case Front:
					localPose.x+=ldelta;
					move=Back;
					break;
				case Back:
					localPose.x-=ldelta;
					move=Left;
					break;
				case Left:
					localPose.y-=ldelta;
					move=Right;
					break;
				case Right:
					localPose.y+=ldelta;
					move=TurnLeft;
					break;
				case TurnLeft:
					localPose.theta+=adelta;
					move=TurnRight;
					break;
				case TurnRight:
					localPose.theta-=adelta;
					move=Done;
					break;
				default:;
			}
			double localScore, localLikelihood;

			double odo_gain=1;
			if (m_angularOdometryReliability>0.){
				double dth=init.theta-localPose.theta; 	dth=atan2(sin(dth), cos(dth)); 	dth*=dth;
				odo_gain*=exp(-m_angularOdometryReliability*dth);
			}
			if (m_linearOdometryReliability>0.){
				double dx=init.x-localPose.x;
				double dy=init.y-localPose.y;
				double drho=dx*dx+dy*dy;
				odo_gain*=exp(-m_linearOdometryReliability*drho);
			}
			localScore=odo_gain*score(map, localPose, readings);
			//update the score
			count++;
			matched=likelihoodAndScore(localScore, localLikelihood, map, localPose, readings);
			if (localScore>currentScore){
				currentScore=localScore;
				bestLocalPose=localPose;
			}
			sm.score=localScore;
			sm.likelihood=localLikelihood;//+log(odo_gain);
			sm.pose=localPose;
			moveList.push_back(sm);
			//update the move list
		} while(move!=Done);
		currentPose=bestLocalPose;
		//cout << __PRETTY_FUNCTION__ << "currentScore=" << currentScore<< endl;
		//here we look for the best move;
	}while (currentScore>bestScore || refinement<m_optRecursiveIterations);
	//cout << __PRETTY_FUNCTION__ << "bestScore=" << bestScore<< endl;
	//cout << __PRETTY_FUNCTION__ << "iterations=" << count<< endl;

	//normalize the likelihood
	double lmin=1e9;
	double lmax=-1e9;
	for (ScoredMoveList::const_iterator it=moveList.begin(); it!=moveList.end(); it++){
		lmin=it->likelihood<lmin?it->likelihood:lmin;
		lmax=it->likelihood>lmax?it->likelihood:lmax;
	}
	//cout << "lmin=" << lmin << " lmax=" << lmax<< endl;
	for (ScoredMoveList::iterator it=moveList.begin(); it!=moveList.end(); it++){
		it->likelihood=exp(it->likelihood-lmax);
		//cout << "l=" << it->likelihood << endl;
	}
	//compute the mean
	OrientedPoint mean(0,0,0);
	double lacc=0;
	for (ScoredMoveList::const_iterator it=moveList.begin(); it!=moveList.end(); it++){
		mean=mean+it->pose*it->likelihood;
		lacc+=it->likelihood;
	}
	mean=mean*(1./lacc);
	//OrientedPoint delta=mean-currentPose;
	//cout << "delta.x=" << delta.x << " delta.y=" << delta.y << " delta.theta=" << delta.theta << endl;
	CovarianceMatrix cov={0.,0.,0.,0.,0.,0.};
	for (ScoredMoveList::const_iterator it=moveList.begin(); it!=moveList.end(); it++){
		OrientedPoint delta=it->pose-mean;
		delta.theta=atan2(sin(delta.theta), cos(delta.theta));
		cov.xx+=delta.x*delta.x*it->likelihood;
		cov.yy+=delta.y*delta.y*it->likelihood;
		cov.tt+=delta.theta*delta.theta*it->likelihood;
		cov.xy+=delta.x*delta.y*it->likelihood;
		cov.xt+=delta.x*delta.theta*it->likelihood;
		cov.yt+=delta.y*delta.theta*it->likelihood;
	}
	cov.xx/=lacc, cov.xy/=lacc, cov.xt/=lacc, cov.yy/=lacc, cov.yt/=lacc, cov.tt/=lacc;

	_mean=currentPose;
	_cov=cov;
	return bestScore;
}

void ScanMatcher::setLaserParameters
	(unsigned int beams, double* angles, const OrientedPoint& lpose){
	/*if (m_laserAngles)
		delete [] m_laserAngles;
	*/
	assert(beams<LASER_MAXBEAMS);
	m_laserPose=lpose;
	m_laserBeams=beams;
	//m_laserAngles=new double[beams];
	memcpy(m_laserAngles, angles, sizeof(double)*m_laserBeams);
}


double ScanMatcher::likelihood
	(double& _lmax, OrientedPoint& _mean, CovarianceMatrix& _cov, const ScanMatcherMap& map, const OrientedPoint& p, const double* readings){
	ScoredMoveList moveList;

	for (double xx=-m_llsamplerange; xx<=m_llsamplerange; xx+=m_llsamplestep)
	for (double yy=-m_llsamplerange; yy<=m_llsamplerange; yy+=m_llsamplestep)
	for (double tt=-m_lasamplerange; tt<=m_lasamplerange; tt+=m_lasamplestep){

		OrientedPoint rp=p;
		rp.x+=xx;
		rp.y+=yy;
		rp.theta+=tt;

		ScoredMove sm;
		sm.pose=rp;

		likelihoodAndScore(sm.score, sm.likelihood, map, rp, readings);
		moveList.push_back(sm);
	}

	//OrientedPoint delta=mean-currentPose;
	//cout << "delta.x=" << delta.x << " delta.y=" << delta.y << " delta.theta=" << delta.theta << endl;
	//normalize the likelihood
	double lmax=-1e9;
	double lcum=0;
	for (ScoredMoveList::const_iterator it=moveList.begin(); it!=moveList.end(); it++){
		lmax=it->likelihood>lmax?it->likelihood:lmax;
	}
	for (ScoredMoveList::iterator it=moveList.begin(); it!=moveList.end(); it++){
		//it->likelihood=exp(it->likelihood-lmax);
		lcum+=exp(it->likelihood-lmax);
		it->likelihood=exp(it->likelihood-lmax);
		//cout << "l=" << it->likelihood << endl;
	}

	OrientedPoint mean(0,0,0);
	double s=0,c=0;
	for (ScoredMoveList::const_iterator it=moveList.begin(); it!=moveList.end(); it++){
		mean=mean+it->pose*it->likelihood;
		s+=it->likelihood*sin(it->pose.theta);
		c+=it->likelihood*cos(it->pose.theta);
	}
	mean=mean*(1./lcum);
	s/=lcum;
	c/=lcum;
	mean.theta=atan2(s,c);


	CovarianceMatrix cov={0.,0.,0.,0.,0.,0.};
	for (ScoredMoveList::const_iterator it=moveList.begin(); it!=moveList.end(); it++){
		OrientedPoint delta=it->pose-mean;
		delta.theta=atan2(sin(delta.theta), cos(delta.theta));
		cov.xx+=delta.x*delta.x*it->likelihood;
		cov.yy+=delta.y*delta.y*it->likelihood;
		cov.tt+=delta.theta*delta.theta*it->likelihood;
		cov.xy+=delta.x*delta.y*it->likelihood;
		cov.xt+=delta.x*delta.theta*it->likelihood;
		cov.yt+=delta.y*delta.theta*it->likelihood;
	}
	cov.xx/=lcum, cov.xy/=lcum, cov.xt/=lcum, cov.yy/=lcum, cov.yt/=lcum, cov.tt/=lcum;

	_mean=mean;
	_cov=cov;
	_lmax=lmax;
	return log(lcum)+lmax;
}

double ScanMatcher::likelihood
	(double& _lmax, OrientedPoint& _mean, CovarianceMatrix& _cov, const ScanMatcherMap& map, const OrientedPoint& p,
	Gaussian3& odometry, const double* readings, double gain){
	ScoredMoveList moveList;


	for (double xx=-m_llsamplerange; xx<=m_llsamplerange; xx+=m_llsamplestep)
	for (double yy=-m_llsamplerange; yy<=m_llsamplerange; yy+=m_llsamplestep)
	for (double tt=-m_lasamplerange; tt<=m_lasamplerange; tt+=m_lasamplestep){

		OrientedPoint rp=p;
		rp.x+=xx;
		rp.y+=yy;
		rp.theta+=tt;

		ScoredMove sm;
		sm.pose=rp;

		likelihoodAndScore(sm.score, sm.likelihood, map, rp, readings);
		sm.likelihood+=odometry.eval(rp)/gain;
		assert(!isnan(sm.likelihood));
		moveList.push_back(sm);
	}

	//OrientedPoint delta=mean-currentPose;
	//cout << "delta.x=" << delta.x << " delta.y=" << delta.y << " delta.theta=" << delta.theta << endl;
	//normalize the likelihood
  double lmax=-std::numeric_limits<double>::max();
	double lcum=0;
	for (ScoredMoveList::const_iterator it=moveList.begin(); it!=moveList.end(); it++){
		lmax=it->likelihood>lmax?it->likelihood:lmax;
	}
	for (ScoredMoveList::iterator it=moveList.begin(); it!=moveList.end(); it++){
		//it->likelihood=exp(it->likelihood-lmax);
		lcum+=exp(it->likelihood-lmax);
		it->likelihood=exp(it->likelihood-lmax);
		//cout << "l=" << it->likelihood << endl;
	}

	OrientedPoint mean(0,0,0);
	double s=0,c=0;
	for (ScoredMoveList::const_iterator it=moveList.begin(); it!=moveList.end(); it++){
		mean=mean+it->pose*it->likelihood;
		s+=it->likelihood*sin(it->pose.theta);
		c+=it->likelihood*cos(it->pose.theta);
	}
	mean=mean*(1./lcum);
	s/=lcum;
	c/=lcum;
	mean.theta=atan2(s,c);


	CovarianceMatrix cov={0.,0.,0.,0.,0.,0.};
	for (ScoredMoveList::const_iterator it=moveList.begin(); it!=moveList.end(); it++){
		OrientedPoint delta=it->pose-mean;
		delta.theta=atan2(sin(delta.theta), cos(delta.theta));
		cov.xx+=delta.x*delta.x*it->likelihood;
		cov.yy+=delta.y*delta.y*it->likelihood;
		cov.tt+=delta.theta*delta.theta*it->likelihood;
		cov.xy+=delta.x*delta.y*it->likelihood;
		cov.xt+=delta.x*delta.theta*it->likelihood;
		cov.yt+=delta.y*delta.theta*it->likelihood;
	}
	cov.xx/=lcum, cov.xy/=lcum, cov.xt/=lcum, cov.yy/=lcum, cov.yt/=lcum, cov.tt/=lcum;

	_mean=mean;
	_cov=cov;
	_lmax=lmax;
	double v=log(lcum)+lmax;
	assert(!isnan(v));
	return v;
}

void ScanMatcher::setMatchingParameters
	(double urange, double range, double sigma, int kernsize, double lopt, double aopt, int iterations,  double likelihoodSigma, unsigned int likelihoodSkip){
	m_usableRange=urange;
	m_laserMaxRange=range;
	m_kernelSize=kernsize;
	m_optLinearDelta=lopt;
	m_optAngularDelta=aopt;
	m_optRecursiveIterations=iterations;
	m_gaussianSigma=sigma;
	m_likelihoodSigma=likelihoodSigma;
	m_likelihoodSkip=likelihoodSkip;
}

};
