#include <errorekf.h>
#include <robot.h>
#include <eigen3/Eigen/Dense>
#include <geometry_msgs/TwistStamped.h>
#include <geometry_msgs/PointStamped.h>
#include <geometry_msgs/PoseStamped.h>
#include <geometry_msgs/WrenchStamped.h>
#include <std_msgs/String.h>
#include <std_msgs/Int32.h>
#include <std_msgs/Bool.h>
#include <mutex>          
#include <thread>        
#include <sensor_msgs/JointState.h>
#include <sensor_msgs/Imu.h>
#include <nav_msgs/Odometry.h>

using namespace Eigen;
using namespace std;

class tracking{
private:
	//ROS Standard Variables
	ros::NodeHandle n;
	//ROS Publishers
	ros::Publisher rel_LLegPose_pub, rel_RLegPose_pub, rel_CoMPose_pub, RLegWrench_pub, LLegWrench_pub, ground_truth_com_pub, ground_truth_odom_pub, legOdom_pub, 
	CoMLegOdom_pub, CoMOdom_pub, COP_pub, baseOdom_pub, RLegOdom_pub, LLegOdom_pub, baseIMU_pub,  ExternalWrench_pub, SupportPose_pub, joint_pub, comp_odom0_pub, SupportLegId_pub;
	//ROS Messages
	nav_msgs::Odometry ground_truth_odom_msg, ground_truth_com_odom_msg, ground_truth_odom_pub_msg, comp_odom0_msg, odom_msg, odom_msg_, ground_truth_odom_msg_;
	//ROS Subscribers
	ros::Subscriber imu_sub, joint_state_sub, lfsr_sub, rfsr_sub, odom_sub, ground_truth_odom_sub, ground_truth_com_sub, compodom0_sub;
	//Buffers for topic data
	Queue<sensor_msgs::JointState> joint_data;
	Queue<sensor_msgs::Imu> base_imu_data;
	Queue<geometry_msgs::WrenchStamped> LLeg_FT_data, RLeg_FT_data;
	std::thread output_thread, filtering_thread;
	std::mutex output_lock;


	Eigen::Affine3d Twv0;
	Eigen::VectorXd joint_state_pos,joint_state_vel;

	std::unique_ptr<Vector3d> wbb, abb, vwb;
	Eigen::Vector3d omegabl, omegabr, vbl, vbr, vbln, vbrn, omegawb, vwl, vwr, omegawr, omegawl, p_FT_LL, p_FT_RL;
	Eigen::Matrix3d JLQnJLt, JRQnJRt;
	Affine3d Twl, Twr, Tbl, Tbr;
	COM_EKF::robotDyn* rd;
    bool useMahony;
	COM_EKF::Madgwick* mw;
    COM_EKF::Mahony* mh;
    double Kp, Ki;
	COM_EKF::deadReckoning* dr;
	Vector3d bias_a,bias_g;
	int imuCalibrationCycles,maxImuCalibrationCycles;
	bool calibrateIMU, computeJointVelocity, data_inc;
   	std::map<std::string, double> joint_state_pos_map, joint_state_vel_map;

	double Tau0, Tau1, VelocityThres;
	double  freq, joint_freq, ft_freq;
	bool odom_inc, check_no_motion;
	bool firstOdom, firstGyrodot, firstJointStates, firstUpdate;
	bool no_motion_indicator, outlier, odom_divergence;
	int number_of_joints, outlier_count;
	bool useGyroLPF;
	int  maWindow;
	int medianWindow;
	int no_motion_it, no_motion_it_threshold;
	double no_motion_threshold;
	std::unique_ptr<Eigen::Quaterniond> q_update;
	Quaterniond  q_update_, q_leg_update, q_now, q_prev;
	std::unique_ptr<Vector3d> pos_update ;
	Vector3d  pos_update_, pos_leg_update, temp, gt_odom;
	Affine3d T_B_A, T_B_G, T_B_P,  , T_FT_LL, T_B_GT;
	Quaterniond  q_B_P, q_B_GT, tempq, qoffsetGTCoM, tempq_, gt_odomq; 
	bool useCoMEKF, useLegOdom, firstGT,firstGTCoM, useOutlierDetection;
    bool debug_mode;


	//Madgwick gain
	double beta;
	// Helper
	bool is_connected_, ground_truth, support_idx_provided;

	Matrix3d Rwb;
	Quaterniond qbs, qbl, qbr, qwb, qwb_, qws, qwl, qwr;
	string base_link_frame, support_foot_frame, lfoot_frame, rfoot_frame;
	


	double mass;
	IMUEKF* imuEKF;
	bool useIMUEKF;
	CoMEKF* nipmEKF;
	butterworthLPF** gyroLPF;
	MovingAverageFilter** gyroMAF;
	//Cuttoff Freqs for LPF
	double gyro_fx, gyro_fy, gyro_fz;
	Vector3d COP_fsr, GRF_fsr, CoM_enc, Gyrodot, Gyro_, CoM_leg_odom;
	double bias_fx,bias_fy,bias_fz;
	JointDF** JointVF;
	double jointFreq,joint_cutoff_freq;
	Mediator *lmdf, *rmdf;

	//WindowMedian<double> *llmdf, *rrmdf;
	string support_leg;

	COM_EKF::ContactDetection* cd;
    Vector3d coplw, coprw;
    double weightl, weightr;


	bool useGEM, ContactDetectionWithCOP, ContactDetectionWithKinematics;
	double foot_polygon_xmin, foot_polygon_xmax, foot_polygon_ymin, foot_polygon_ymax;
	double lforce_sigma, rforce_sigma, lcop_sigma, rcop_sigma, lvnorm_sigma, rvnorm_sigma, probabilisticContactThreshold;
	Vector3d LLegGRF, RLegGRF, LLegGRT, RLegGRT, offsetGT,offsetGTCoM;
  	Vector3d copl, copr;
	Affine3d Tws, Twb, Twb_; //From support s to world frame;
	Affine3d Tbs, Tsb, Tssw, Tbsw;
	Vector3d no_motion_residual;
	/****/
	bool  kinematicsInitialized, firstContact;
	Vector3d LLegForceFilt, RLegForceFilt;
	double LegHighThres, LegLowThres, LosingContact, StrikingContact;
	double bias_ax, bias_ay, bias_az, bias_gx, bias_gy, bias_gz;
	double g, I_xx, I_yy, I_zz;
	double joint_noise_density;

	bool comp_with, comp_odom0_inc, firstCO;
	std::string comp_with_odom0_topic;
	Vector3d offsetCO;
	Quaterniond qoffsetCO;
	void subscribeToCompOdom();
	void compodom0Cb(const nav_msgs::Odometry::ConstPtr& msg);
	/** Real odometry Data **/
     string lfsr_topic,rfsr_topic; 
	 string imu_topic;
	 string joint_state_topic;
	 string odom_topic;
	 string ground_truth_odom_topic, ground_truth_com_topic;
     string modelname;
	//Odometry, from supportleg to inertial, transformation from support leg to other leg
     void subscribeToIMU();
	 void subscribeToFSR();
	 void subscribeToJointState();
	 void subscribeToOdom();
	 void subscribeToGroundTruth();
	 void subscribeToGroundTruthCoM();
	 void ground_truth_comCb(const nav_msgs::Odometry::ConstPtr& msg);
	 void ground_truth_odomCb(const nav_msgs::Odometry::ConstPtr& msg);
	 void imuCb(const sensor_msgs::Imu::ConstPtr& msg);
	 void joint_stateCb(const sensor_msgs::JointState::ConstPtr& msg);
	 void odomCb(const nav_msgs::Odometry::ConstPtr& msg);
	 void lfsrCb(const geometry_msgs::WrenchStamped::ConstPtr& msg);
	 void rfsrCb(const geometry_msgs::WrenchStamped::ConstPtr& msg);
	 void computeGlobalCOP(Affine3d Tis_, Affine3d Tssprime_);
	 void filterGyrodot();
	//private methods
	void init();
	void estimateWithCoMEKF();
	void estimateWithIMUEKF();
	void computeKinTFs();
	//publish functions
	void publishGRF();
	void publishJointEstimates();
	void publishCoMEstimates();
	void deAllocate();
	void publishLegEstimates();
	void publishSupportEstimates();
	void publishBodyEstimates();
	void publishContact();
	void publishCOP();
	// Advertise to ROS Topics
	void advertise();
	void subscribe();
	void outputPublishThread();
	void filteringThread();
	void joints(const sensor_msgs::JointState &msg);
	void baseIMU(const sensor_msgs::Imu &msg);
	void LLeg_FT(const geometry_msgs::WrenchStamped &msg);
	void RLeg_FT(const geometry_msgs::WrenchStamped &msg);
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	// Constructor/Destructor
	humanoid_ekf();
	~humanoid_ekf();
	bool connect(const ros::NodeHandle nh);
	void disconnect();
	// Parameter Server
	void loadparams();
	void loadIMUEKFparams();
	void loadCoMEKFparams();
	void loadJointKFparams();
	void run();
	bool connected();
};

