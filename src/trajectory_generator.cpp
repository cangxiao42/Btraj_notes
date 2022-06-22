#include "trajectory_generator.h"
using namespace std;    
using namespace Eigen;

static void MSKAPI printstr(void *handle, MSKCONST char str[])
{
  printf("%s",str);
}

int TrajectoryGenerator::BezierPloyCoeffGeneration(
            const vector<Cube> &corridor,
            const MatrixXd &MQM,
            const MatrixXd &pos,
            const MatrixXd &vel,
            const MatrixXd &acc,
            const double maxVel,
            const double maxAcc,
            const int traj_order,
            const double minimize_order,
            const double margin,
            const bool & isLimitVel,
            const bool & isLimitAcc,
            double & obj,
            MatrixXd & PolyCoeff)  // define the order to which we minimize.   1 -- velocity, 2 -- acceleration, 3 -- jerk, 4 -- snap  
{   
#define ENFORCE_VEL  isLimitVel // whether or not adding extra constraints for ensuring the velocity feasibility
#define ENFORCE_ACC  isLimitAcc // whether or not adding extra constraints for ensuring the acceleration feasibility

    double initScale = corridor.front().t;
    double lstScale  = corridor.back().t;
    int segment_num  = corridor.size();

    int n_poly = traj_order + 1; // 多项式系数的个数等于阶数加一
    int s1d1CtrlP_num = n_poly; // 一段轨迹的一维坐标的控制点的数量
    int s1CtrlP_num   = 3 * s1d1CtrlP_num; // 一段轨迹中三个维度的控制点总数

    int equ_con_s_num = 3 * 3; // p, v, a in x, y, z axis at the start point，起点约束
    int equ_con_e_num = 3 * 3; // p, v, a in x, y, z axis at the end point，终点约束
    int equ_con_continuity_num = 3 * 3 * (segment_num - 1); // 三个方向的p，v，a连续性约束
    int equ_con_num   = equ_con_s_num + equ_con_e_num + equ_con_continuity_num; // p, v, a in x, y, z axis in each segment's joint position，总的等式约束是3n+3个，之所以不是4n+2个，是因为没有对路径点进行位置约束，所以少了n-1个约束
    
    /**
     * @param vel_con_num 在所有的控制点上速度都要受到速度约束，这里是不等式约束，之所以是(3 * traj_order * segment_num)是因为，三个维度，segment_num段轨迹，
     * 贝赛尔曲线有（traj_order+1）个控制点,控制轨迹的位置，而贝赛尔曲线的特性之一是，曲线的导数还是贝赛尔曲线
     * 而且，导数的控制点等于n（控制点（i） - 控制点（i-1））所以一阶导的控制点是（traj_order+1）- 1个，也就是traj_order个控制点
     * 同理，二阶导的控制点有traj_order-1个
     */
    int vel_con_num = 3 *  traj_order * segment_num;  
    int acc_con_num = 3 * (traj_order - 1) * segment_num; 

    if( !ENFORCE_VEL )
        vel_con_num = 0;

    if( !ENFORCE_ACC )
        acc_con_num = 0;

    int high_order_con_num = vel_con_num + acc_con_num; 
    //int high_order_con_num = 0; //3 * traj_order * segment_num;

    int con_num   = equ_con_num + high_order_con_num; // 所有轨迹，所有维度的约束的个数，包含等式和不等式约束
    int ctrlP_num = segment_num * s1CtrlP_num; // 所有轨迹，所有维度的控制点的个数

    double x_var[ctrlP_num]; // 存放控制点的最终解
    double primalobj;

    MSKrescodee  r; 
    /**
     * @brief 构建所有关于constraint的约束条件，包含等式约束和不等式约束，等式约束包含起点约束，终点约束，连续性约束，不等式约束包含速度和加速度的约束
     * @param con_bdk 中，MSKboundkeye是约束的键值，表示约束的类型，具体可参考Mosek官方手册：https://docs.mosek.com/8.1/capi/conventions.html#doc-optimizer-cmo-rmo-matrix
     * 
     */
    vector< pair<MSKboundkeye, pair<double, double> > > con_bdk; 
    
    if(ENFORCE_VEL)
    {
        /***  Stack the bounding value for the linear inequality for the velocity constraints  ***/
        for(int i = 0; i < vel_con_num; i++)
        {
            pair<MSKboundkeye, pair<double, double> > cb_ie = make_pair( MSK_BK_RA, make_pair( - maxVel,  + maxVel) );
            con_bdk.push_back(cb_ie);   
        }
    }

    if(ENFORCE_ACC)
    {
        /***  Stack the bounding value for the linear inequality for the acceleration constraints  ***/
        for(int i = 0; i < acc_con_num; i++)
        {
            pair<MSKboundkeye, pair<double, double> > cb_ie = make_pair( MSK_BK_RA, make_pair( - maxAcc,  maxAcc) ); 
            con_bdk.push_back(cb_ie);   
        }
    }

    //ROS_WARN("[Bezier Trajectory] equality bound %d", equ_con_num);
    // 这里主要是约束起点终点的pva，Ax=b中，起点终点的pva是硬约束，b中需要定值，而其他的约束的对应的b为0即可
    for(int i = 0; i < equ_con_num; i ++ ){ 
        double beq_i;
        if(i < 3)                    beq_i = pos(0, i); 
        else if (i >= 3  && i < 6  ) beq_i = vel(0, i - 3); 
        else if (i >= 6  && i < 9  ) beq_i = acc(0, i - 6);
        else if (i >= 9  && i < 12 ) beq_i = pos(1, i - 9 );
        else if (i >= 12 && i < 15 ) beq_i = vel(1, i - 12);
        else if (i >= 15 && i < 18 ) beq_i = acc(1, i - 15);
        else beq_i = 0.0;

        pair<MSKboundkeye, pair<double, double> > cb_eq = make_pair( MSK_BK_FX, make_pair( beq_i, beq_i ) ); // # cb_eq means: constriants boundary of equality constrain
        con_bdk.push_back(cb_eq);
    }

    /* ## define a container for control points' boundary and boundkey ## */ 
    /* ## dataType in one tuple is : boundary type, lower bound, upper bound ## */
    /**
     * @brief 此处计算关于控制点自身的约束，也就是变量的约束，不等式约束包含控制点的位置约束
     * @param margin 应该是设定控制点的边界，因为普通的box是顶在障碍物上的，这里是要与box的边界保持一定距离
     * @param scale_k 表示放缩的参数，用来对时间进行归一话，具体可参照论文，论文中有详细描述
     * 这里对控制点不需要对轨迹连接处的控制点单独约束是因为，控制点普通方法约束就可以将其约束在两个box的重叠区域
     * 例如：           
     *             ____________
     *            |            |box1
     *            |            |
     *            |       _____|_______     
     *            |      |ctrl |       |     
     *            |______|_____|       |
     *                   |             |
     *                   |             |
     *                   |_____ _______|
     *                          box2
     * ctrl表示一个控制点，ctrl即约束在box1内，又约束在box2内，就相当于ctrl被约束在了两个box的重叠区域，所以不需要单独计算
     * 两个box的重叠区域作为约束了
     */ 

    vector< pair<MSKboundkeye, pair<double, double> > > var_bdk; 

    for(int k = 0; k < segment_num; k++)
    {   
        Cube cube_     = corridor[k];
        double scale_k = cube_.t;

        for(int i = 0; i < 3; i++ )
        {   
            for(int j = 0; j < n_poly; j ++ )
            {   
                pair<MSKboundkeye, pair<double, double> > vb_x;

                double lo_bound, up_bound;
                if(k > 0)
                {
                    lo_bound = (cube_.box[i].first  + margin) / scale_k;
                    up_bound = (cube_.box[i].second - margin) / scale_k;
                }
                else
                {
                    lo_bound = (cube_.box[i].first)  / scale_k;
                    up_bound = (cube_.box[i].second) / scale_k;
                }

                vb_x  = make_pair( MSK_BK_RA, make_pair( lo_bound, up_bound ) ); // # vb_x means: varialbles boundary of unknowns x (Polynomial coeff)

                var_bdk.push_back(vb_x);
            }
        } 
    }

    /*****************************************建立Mosek的相关配置和环境*****************************************************/
    MSKint32t  j,i; 
    MSKenv_t   env; 
    MSKtask_t  task; 
    // Create the mosek environment. 
    r = MSK_makeenv( &env, NULL ); 
  
    // Create the optimization task. 
    r = MSK_maketask(env,con_num, ctrlP_num, &task); // con_num是所有constraint的数量，ctrlP_num是所有控制点的数量

// Parameters used in the optimizer
//######################################################################
    //MSK_putintparam (task, MSK_IPAR_OPTIMIZER , MSK_OPTIMIZER_INTPNT );
    MSK_putintparam (task, MSK_IPAR_NUM_THREADS, 1);
    MSK_putdouparam (task, MSK_DPAR_CHECK_CONVEXITY_REL_TOL, 1e-2);
    MSK_putdouparam (task, MSK_DPAR_INTPNT_TOL_DFEAS,  1e-4);
    MSK_putdouparam (task, MSK_DPAR_INTPNT_TOL_PFEAS,  1e-4);
    MSK_putdouparam (task, MSK_DPAR_INTPNT_TOL_INFEAS, 1e-4);
    //MSK_putdouparam (task, MSK_DPAR_INTPNT_TOL_REL_GAP, 5e-2 );
//######################################################################
    
    //r = MSK_linkfunctotaskstream(task,MSK_STREAM_LOG,NULL,printstr); 
    // Append empty constraints. 
     //The constraints will initially have no bounds. 
     /*********************************************开始构建Ax=b中的A************************************************/
    if ( r == MSK_RES_OK ) 
      r = MSK_appendcons(task,con_num);  

    // Append optimizing variables. The variables will initially be fixed at zero (x=0). 
    if ( r == MSK_RES_OK ) 
      r = MSK_appendvars(task,ctrlP_num); 

    //ROS_WARN("set variables boundary");
    for(j = 0; j<ctrlP_num && r == MSK_RES_OK; ++j){ 
        if (r == MSK_RES_OK) 
            r = MSK_putvarbound(task, 
                                j,                            // Index of variable. 
                                var_bdk[j].first,             // Bound key.
                                var_bdk[j].second.first,      // Numerical value of lower bound.
                                var_bdk[j].second.second );   // Numerical value of upper bound.      
    } 
    
    // Set the bounds on constraints. 
    //   for i=1, ...,con_num : blc[i] <= constraint i <= buc[i] 
    for( i = 0; i < con_num && r == MSK_RES_OK; i++ ) {
        r = MSK_putconbound(task, 
                            i,                            // Index of constraint. 
                            con_bdk[i].first,             // Bound key.
                            con_bdk[i].second.first,      // Numerical value of lower bound.
                            con_bdk[i].second.second );   // Numerical value of upper bound. 
    }

    //ROS_WARN("[Bezier Trajectory] Start stacking the Linear Matrix A, inequality part");
    int row_idx = 0;
    // The velocity constraints
     /**
     * @brief 关于 @param MSK_putarow 的用法: https://docs.mosek.com/8.1/capi/alphabetic-functionalities.html#mosek.task.putarow
     * 这里做简单介绍
     * 此函数首先将矩阵的第row_idx行全部置为0
     * @param nzi 接下来需要将第row_idx行某些列变成非0数，nzi表示非0数的数量
     * @param asub 存放需要转变的非0数的列号
     * @param aval 存放对应位置的新值
     * 比如，要将第3行第4列变为111，第五列变为222
     * 那么 row_idx = 3 
     * nzi = 2
     * asub = [4, 5]
     * aval = [111, 222]
     * 这么表示主要是考虑稀疏矩阵的表示形式
     */
    if(ENFORCE_VEL)
    {   
        for(int k = 0; k < segment_num ; k ++ )
        {   
            for(int i = 0; i < 3; i++)
            {  // for x, y, z loop
                for(int p = 0; p < traj_order; p++)
                {
                    int nzi = 2;
                    MSKint32t asub[nzi];
                    double aval[nzi];

                    aval[0] = -1.0 * traj_order; // traj_order*(控制点[1] - 控制点[i-1])表示有贝赛尔曲线得到的对应的速度
                    aval[1] =  1.0 * traj_order;

                    // 对应的控制点的位置，由于控制点就是优化变量，而且以向量形式保存，所以计算位置
                    // 还原过来就是 v = 控制点[asub[1]]*aval[1] + 控制点[asub[0]]*aval[0]，用+号是因为aval[0]已经是负的了
                    asub[0] = k * s1CtrlP_num + i * s1d1CtrlP_num + p;    
                    asub[1] = k * s1CtrlP_num + i * s1d1CtrlP_num + p + 1;    

                    r = MSK_putarow(task, row_idx, nzi, asub, aval);    
                    row_idx ++;
                }
            }
        }
    }

    // The acceleration constraints
    if(ENFORCE_ACC)
    {
        for(int k = 0; k < segment_num ; k ++ )
        {
            for(int i = 0; i < 3; i++)
            { 
                for(int p = 0; p < traj_order - 1; p++)
                {    
                    int nzi = 3;
                    MSKint32t asub[nzi];
                    double aval[nzi];
                    // 此处与上面同理
                    aval[0] =  1.0 * traj_order * (traj_order - 1) / corridor[k].t;
                    aval[1] = -2.0 * traj_order * (traj_order - 1) / corridor[k].t;
                    aval[2] =  1.0 * traj_order * (traj_order - 1) / corridor[k].t;
                    asub[0] = k * s1CtrlP_num + i * s1d1CtrlP_num + p;    
                    asub[1] = k * s1CtrlP_num + i * s1d1CtrlP_num + p + 1;    
                    asub[2] = k * s1CtrlP_num + i * s1d1CtrlP_num + p + 2;    
                    
                    r = MSK_putarow(task, row_idx, nzi, asub, aval);    
                    row_idx ++;
                }
            }
        }
    }
    // 下面这些设置与上面同理，不明白的地方可以看论文
    /*   Start position  */
    {
        // position :
        for(int i = 0; i < 3; i++)
        {  // loop for x, y, z       
            int nzi = 1;
            MSKint32t asub[nzi];
            double aval[nzi];
            aval[0] = 1.0 * initScale;
            asub[0] = i * s1d1CtrlP_num;
            r = MSK_putarow(task, row_idx, nzi, asub, aval);    
            row_idx ++;
        }
        // velocity :
        for(int i = 0; i < 3; i++)
        {  // loop for x, y, z       
            int nzi = 2;
            MSKint32t asub[nzi];
            double aval[nzi];
            aval[0] = - 1.0 * traj_order;
            aval[1] =   1.0 * traj_order;
            asub[0] = i * s1d1CtrlP_num;
            asub[1] = i * s1d1CtrlP_num + 1;
            r = MSK_putarow(task, row_idx, nzi, asub, aval);   
            row_idx ++;
        }
        // acceleration : 
        for(int i = 0; i < 3; i++)
        {  // loop for x, y, z       
            int nzi = 3;
            MSKint32t asub[nzi];
            double aval[nzi];
            aval[0] =   1.0 * traj_order * (traj_order - 1) / initScale;
            aval[1] = - 2.0 * traj_order * (traj_order - 1) / initScale;
            aval[2] =   1.0 * traj_order * (traj_order - 1) / initScale;
            asub[0] = i * s1d1CtrlP_num;
            asub[1] = i * s1d1CtrlP_num + 1;
            asub[2] = i * s1d1CtrlP_num + 2;
            r = MSK_putarow(task, row_idx, nzi, asub, aval);    
            row_idx ++;
        }
    }      

    /*   End position  */
    //ROS_WARN(" end position");
    {   
        // position :
        for(int i = 0; i < 3; i++)
        {  // loop for x, y, z       
            int nzi = 1;
            MSKint32t asub[nzi];
            double aval[nzi];
            asub[0] = ctrlP_num - 1 - (2 - i) * s1d1CtrlP_num;
            aval[0] = 1.0 * lstScale;
            r = MSK_putarow(task, row_idx, nzi, asub, aval);    
            row_idx ++;
        }
        // velocity :
        for(int i = 0; i < 3; i++)
        { 
            int nzi = 2;
            MSKint32t asub[nzi];
            double aval[nzi];
            asub[0] = ctrlP_num - 1 - (2 - i) * s1d1CtrlP_num - 1;
            asub[1] = ctrlP_num - 1 - (2 - i) * s1d1CtrlP_num;
            aval[0] = - 1.0;
            aval[1] =   1.0;
            r = MSK_putarow(task, row_idx, nzi, asub, aval);    
            row_idx ++;
        }
        // acceleration : 
        for(int i = 0; i < 3; i++)
        { 
            int nzi = 3;
            MSKint32t asub[nzi];
            double aval[nzi];
            asub[0] = ctrlP_num - 1 - (2 - i) * s1d1CtrlP_num - 2;
            asub[1] = ctrlP_num - 1 - (2 - i) * s1d1CtrlP_num - 1;
            asub[2] = ctrlP_num - 1 - (2 - i) * s1d1CtrlP_num;
            aval[0] =   1.0 / lstScale;
            aval[1] = - 2.0 / lstScale;
            aval[2] =   1.0 / lstScale;
            r = MSK_putarow(task, row_idx, nzi, asub, aval);    
            row_idx ++;
        }
    }

    /*   joint points  */
    //ROS_WARN(" joint position");
    // 这里构建各段轨迹的连续性约束，并在A中表示
    // 比如，第一段轨迹末尾和第二段轨迹的起点应该保持位置，速度，加速度等的连续
    {
        int sub_shift = 0;
        double val0, val1;
        for(int k = 0; k < (segment_num - 1); k ++ )
        {   // 这里前后两段轨迹放缩程度不同
            double scale_k = corridor[k].t;
            double scale_n = corridor[k+1].t;
            // position :
            val0 = scale_k;
            val1 = scale_n;
            // 位置连续
            for(int i = 0; i < 3; i++)
            {  // loop for x, y, z
                int nzi = 2;
                MSKint32t asub[nzi];
                double aval[nzi];

                // This segment's last control point
                aval[0] = 1.0 * val0;
                asub[0] = sub_shift + (i+1) * s1d1CtrlP_num - 1;

                // Next segment's first control point
                aval[1] = -1.0 * val1;
                asub[1] = sub_shift + s1CtrlP_num + i * s1d1CtrlP_num;
                r = MSK_putarow(task, row_idx, nzi, asub, aval);    
                row_idx ++;
            }
            
            // 速度连续
            for(int i = 0; i < 3; i++)
            {  
                int nzi = 4;
                MSKint32t asub[nzi];
                double aval[nzi];
                
                // This segment's last velocity control point
                aval[0] = -1.0;
                aval[1] =  1.0;
                asub[0] = sub_shift + (i+1) * s1d1CtrlP_num - 2;    
                asub[1] = sub_shift + (i+1) * s1d1CtrlP_num - 1;   
                // Next segment's first velocity control point
                aval[2] =  1.0;
                aval[3] = -1.0;

                asub[2] = sub_shift + s1CtrlP_num + i * s1d1CtrlP_num;    
                asub[3] = sub_shift + s1CtrlP_num + i * s1d1CtrlP_num + 1;

                r = MSK_putarow(task, row_idx, nzi, asub, aval);    
                row_idx ++;
            }
            // acceleration :
            val0 = 1.0 / scale_k;
            val1 = 1.0 / scale_n;
            // 加速度连续
            for(int i = 0; i < 3; i++)
            {  
                int nzi = 6;
                MSKint32t asub[nzi];
                double aval[nzi];
                
                // This segment's last velocity control point
                aval[0] =  1.0  * val0;
                aval[1] = -2.0  * val0;
                aval[2] =  1.0  * val0;
                asub[0] = sub_shift + (i+1) * s1d1CtrlP_num - 3;    
                asub[1] = sub_shift + (i+1) * s1d1CtrlP_num - 2;   
                asub[2] = sub_shift + (i+1) * s1d1CtrlP_num - 1;   
                // Next segment's first velocity control point
                aval[3] =  -1.0  * val1;
                aval[4] =   2.0  * val1;
                aval[5] =  -1.0  * val1;
                asub[3] = sub_shift + s1CtrlP_num + i * s1d1CtrlP_num;    
                asub[4] = sub_shift + s1CtrlP_num + i * s1d1CtrlP_num + 1;
                asub[5] = sub_shift + s1CtrlP_num + i * s1d1CtrlP_num + 2;

                r = MSK_putarow(task, row_idx, nzi, asub, aval);    
                row_idx ++;
            }

            sub_shift += s1CtrlP_num;
        }
    }

    //ROS_WARN("[Bezier Trajectory] Start stacking the objective");
    
    /*****************************构造优化目标函数的矩阵M‘QM****************************/
    // Q是minimum-jerk或minimum-snap得到的二次形，而M是为了将普通多项式转变为贝赛尔多项式构造的矩阵
    // 原本的优化函数形式：p’Qp，p表示原本的多项式系数，没有实际的物理含义
    // 令p = Mu， u表示贝赛尔曲线的控制点，这样u就有了物理意义，优化函数转变为u‘M’QMu
    int min_order_l = floor(minimize_order);
    int min_order_u = ceil (minimize_order);

    /**
     * @brief Q矩阵是对称矩阵，所以仅以下三角的形式表示即可
     * @param NUMQNZ 计算表示Q矩阵所必需的元素的个数
     * 例如：
     * 七次多项式，有八个系数，所以Q是8×8的矩阵，由于其对称，所以我们取下三角来表示，那么元素的个数是：
     * （（8×8）- 8）/2 + 8
     * 8×8是元素的总个数，-8是减去对角线元素，除2之后得到不包含对角线元素的下三角的元素的个数，加8是加上对角线的元素，花间之后就得到：
     * NUMQ_blk * (NUMQ_blk + 1) / 2
     */
    int NUMQNZ = 0;
    for(int i = 0; i < segment_num; i ++)
    {   
        int NUMQ_blk = (traj_order + 1);                       // default minimize the jerk and minimize_order = 3
        NUMQNZ      += 3 * NUMQ_blk * (NUMQ_blk + 1) / 2; 
    }
    MSKint32t  qsubi[NUMQNZ], qsubj[NUMQNZ];
    double     qval[NUMQNZ];
    
    {    
        int sub_shift = 0;
        int idx = 0;
        for(int k = 0; k < segment_num; k ++)
        {
            double scale_k = corridor[k].t;
            for(int p = 0; p < 3; p ++ )
                for( int i = 0; i < s1d1CtrlP_num; i ++ )
                    for( int j = 0; j < s1d1CtrlP_num; j ++ )
                        if( i >= j )
                        {
                            // 将矩阵转换为3个数组, qsubi: 每个元素的行数, qsubj: 每个元素的列数, qval: 每个元素的值
                            qsubi[idx] = sub_shift + p * s1d1CtrlP_num + i;   
                            qsubj[idx] = sub_shift + p * s1d1CtrlP_num + j;  
                            //qval[idx]  = MQM(i, j) /(double)pow(scale_k, 3);
                            if(min_order_l == min_order_u)
                                qval[idx]  = MQM(i, j) /(double)pow(scale_k, 2 * min_order_u - 3);
                            else
                                qval[idx] = ( (minimize_order - min_order_l) / (double)pow(scale_k, 2 * min_order_u - 3)
                                            + (min_order_u - minimize_order) / (double)pow(scale_k, 2 * min_order_l - 3) ) * MQM(i, j);
                            idx ++ ;
                        }

            sub_shift += s1CtrlP_num;
        }
    }
         
    ros::Time time_end1 = ros::Time::now();
    /******************************Mosek求解环节**********************************/
    if ( r== MSK_RES_OK )
         r = MSK_putqobj(task,NUMQNZ,qsubi,qsubj,qval); 
    
    if ( r==MSK_RES_OK ) 
         r = MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE);
    
    //ros::Time time_opt = ros::Time::now();
    bool solve_ok = false;
    if ( r==MSK_RES_OK ) 
      { 
        //ROS_WARN("Prepare to solve the problem ");   
        MSKrescodee trmcode; 
        r = MSK_optimizetrm(task,&trmcode); 
        MSK_solutionsummary (task,MSK_STREAM_LOG); 
          
        if ( r==MSK_RES_OK ) 
        { 
          MSKsolstae solsta; 
          MSK_getsolsta (task,MSK_SOL_ITR,&solsta); 
           
          switch(solsta) 
          { 
            case MSK_SOL_STA_OPTIMAL:    
            case MSK_SOL_STA_NEAR_OPTIMAL: 
              
            
            r = MSK_getxx(task, 
                          MSK_SOL_ITR,    // Request the interior solution.  
                          x_var); 
            
            r = MSK_getprimalobj(
                task,
                MSK_SOL_ITR,
                &primalobj);

            obj = primalobj;
            solve_ok = true;
            
            break; 
            
            case MSK_SOL_STA_DUAL_INFEAS_CER: 
            case MSK_SOL_STA_PRIM_INFEAS_CER: 
            case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER: 
            case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:   
              printf("Primal or dual infeasibility certificate found.\n"); 
              break; 
               
            case MSK_SOL_STA_UNKNOWN: 
              printf("The status of the solution could not be determined.\n"); 
              //solve_ok = true; // debug
              break; 
            default: 
              printf("Other solution status."); 
              break; 
          } 
        } 
        else 
        { 
          printf("Error while optimizing.\n"); 
        } 
      }
     
      if (r != MSK_RES_OK) 
      { 
        // In case of an error print error code and description. 
        char symname[MSK_MAX_STR_LEN]; 
        char desc[MSK_MAX_STR_LEN]; 
         
        printf("An error occurred while optimizing.\n");      
        MSK_getcodedesc (r, 
                         symname, 
                         desc); 
        printf("Error %s - '%s'\n",symname,desc); 
      } 
    
    MSK_deletetask(&task); 
    MSK_deleteenv(&env); 

    ros::Time time_end2 = ros::Time::now();
    ROS_WARN("time consume in optimize is :");
    cout<<time_end2 - time_end1<<endl;

    if(!solve_ok){
      ROS_WARN("In solver, falied ");
      return -1;
    }
    /*************************输出结果，将控制点存放导PloyCoeff****************************/
    VectorXd d_var(ctrlP_num);
    for(int i = 0; i < ctrlP_num; i++)
        d_var(i) = x_var[i];
    
    PolyCoeff = MatrixXd::Zero(segment_num, 3 *(traj_order + 1) );

    int var_shift = 0;
    for(int i = 0; i < segment_num; i++ )
    {
        for(int j = 0; j < 3 * n_poly; j++)
            PolyCoeff(i , j) = d_var(j + var_shift);

        var_shift += 3 * n_poly;
    }   

    return 1;
}