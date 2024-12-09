#include "kalman_filter.h"
void KalmanFilter_Horizontal_XYZ(float mea_pos,         //位置测量值(m)
								 float mea_vel,         //速度测量值(m/s)
								 float mea_acc,         //加速度测量值(m/s^2)
								 float *optimal_est_pos,//最优估算位置(m)
								 float *optimal_est_vel,//最优估算速度(m/s)
								 uint8_t Axis,          //选择坐标轴
								 float dt)              //dt时间间隔执行一次(s)
{
	uint8_t Label = 0;
	float kalman_A[2][2] = {{1  ,  dt}, //状态空间矩阵
	                        {0  ,  1 }};
	double conv_matrix_pre[2][2] = {{0  ,  0},
									{0  ,  0}}; //先验协方差矩阵
	double kalman_gain[2][2] = {{0  ,  0},//卡尔曼增益矩阵
                              	{0  ,  0}};
	float conv_z = 0;//求卡尔曼增益分母上的逆矩阵时的行列式的值
	float zdelta[2] = {0};
	if(Axis == 'X') Label = 0;
	else if(Axis == 'Y')Label = 1;
	else if(Axis == 'Z')Label = 2;
	//先验估计
	/*
	标准形式：Xk~ = AX~k-1 + BUk-1
	状态空间方程表达式
	pk~ = pk-1 + vk-1 * dt + 0.5dt^2
	vk~ = vk-1 + ak-1 * dt
	状态空间方程矩阵形式
	[p]k   [1   dt][p]k-1    [0.5dt^2]
		 =                 +          [a]k-1
	[v]k   [0    1][v]k-1    [   dt  ]

				[1   dt]    [0.5dt^2]
	其中A=          B=
				[0    1]    [  dt   ]
	*/
	*optimal_est_pos = kalman_A[0][0] * *optimal_est_pos + kalman_A[0][1] * *optimal_est_vel + 0.5f * mea_acc * dt * dt;
	*optimal_est_vel = kalman_A[1][0] * *optimal_est_pos + kalman_A[1][1] * *optimal_est_vel + mea_acc * dt;
	//先验协方差矩阵
	/*
	标准形式：Pk~ = APk-1A^T + Q
	先验误差协方差矩阵 = 状态空间矩阵 * 上一次的协方差矩阵 * 状态空间矩阵的转置 + 过程噪声的协方差矩阵
	*/
	conv_matrix_pre[0][0] =  conv_matrix[Label][0][0] + conv_matrix[Label][1][0] * dt + (conv_matrix[Label][0][1] + conv_matrix[Label][1][1] * dt) * dt + Q[Label][0][0];
	conv_matrix_pre[0][1] =  conv_matrix[Label][0][1] + conv_matrix[Label][1][1] * dt + Q[Label][0][1];
	conv_matrix_pre[1][0] =  conv_matrix[Label][1][0] + conv_matrix[Label][1][1] * dt + Q[Label][1][0];
	conv_matrix_pre[1][1] =  conv_matrix[Label][1][1] + Q[Label][1][1];
	//计算卡尔曼增益(先对分母上的矩阵求逆矩阵)
	/*
										Pk~H^T
	标准形式：Kk = ——————
									HPk~H^T + R
	卡尔曼增益 = （先验误差协方差矩阵Pk~ * 状态转移矩阵H的转置）* （状态转移矩阵H * 先验误差协方差矩阵Pk~ * 状态转移矩阵H的转置 + 测量噪声的协方差矩阵R）^-1
	本实验中，状态转移矩阵H为单位阵E，即公式可简化为

						Pk~
	Kk = ——————
							Pk~ + R
		等价位Pk~ * (Pk~ + R)的逆
	*/
	conv_z = 1.0f / ((conv_matrix_pre[0][0] + R[Label][0][0]) * (conv_matrix_pre[1][1] + R[Label][1][1]) - (conv_matrix_pre[0][1] + R[Label][0][1]) * (conv_matrix_pre[1][0] + R[Label][1][0]));
	kalman_gain[0][0] = conv_z * ( conv_matrix_pre[0][0] * (conv_matrix_pre[1][1] + R[Label][1][1]) - conv_matrix_pre[0][1] * (conv_matrix_pre[1][0] + R[Label][1][0]));
	kalman_gain[0][1] = conv_z * (-conv_matrix_pre[0][0] * R[Label][0][1] + conv_matrix_pre[0][1] * R[Label][0][0]);
	kalman_gain[1][0] = conv_z * ( conv_matrix_pre[1][0] * R[Label][1][1] - conv_matrix_pre[1][1] * R[Label][1][0]);
	kalman_gain[1][1] = conv_z * (-conv_matrix_pre[1][0] * (conv_matrix_pre[0][1] + R[Label][0][1]) + conv_matrix_pre[1][1] * (conv_matrix_pre[0][0] + R[Label][0][0]));	
	//后验估计
	/*
	标准形式：Xk = Xk~ + Kk(Zk - HXk~)
	后验估计（最优估计） = 先验估计 + 卡尔曼增益矩阵 * （测量值 - 状态转移矩阵 * 先验估计）

	[p]k   [p]k~     [Kk00 Kk01][(Zpk-1) - pk~] 
			 =        +                
	[v]k   [v]k~     [Kk10 Kk11][(Zvk-1) - vk~]
	*/
	zdelta[0] = constrainf(mea_pos - *optimal_est_pos,-5,5);
	zdelta[1] = constrainf(mea_vel - *optimal_est_vel,-5,5);
	*optimal_est_pos = *optimal_est_pos + kalman_gain[0][0] * zdelta[0] + kalman_gain[0][1] * zdelta[1];
	*optimal_est_vel = *optimal_est_vel + kalman_gain[1][0] * zdelta[0] + kalman_gain[1][1] * zdelta[1];
	//更新状态协方差矩阵	
	/*
	标准形式:  Pk = (I - Kk*H)*Pk~
	状态协方差矩阵 = （单位阵 - 卡尔曼增益矩阵 * 状态转移矩阵） * 先验协方差矩阵
	*/		
	conv_matrix[Label][0][0] = (1 - kalman_gain[0][0]) * conv_matrix_pre[0][0] - kalman_gain[0][1] * conv_matrix_pre[1][0];
	conv_matrix[Label][0][1] = (1 - kalman_gain[0][0]) * conv_matrix_pre[0][1] - kalman_gain[0][1] * conv_matrix_pre[1][1];
	conv_matrix[Label][1][0] = (1 - kalman_gain[1][1]) * conv_matrix_pre[1][0] - kalman_gain[1][0] * conv_matrix_pre[0][0];
	conv_matrix[Label][1][1] = (1 - kalman_gain[1][1]) * conv_matrix_pre[1][1] - kalman_gain[1][0] * conv_matrix_pre[0][1];
}