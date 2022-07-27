#include "mws_ivp_solver.h"
#include "mws_ls_solver.h"
#include "mws_nls_solver.h"

#include "my_RK45.c"    /*�Զ����㷨ͷ�ļ�*/

void MwsRegisterUserAlgorithm1(void* mdl_data)
{
	/* Register user defined LS and NLS algorithm. ���Ի�������㷨*/
}

void MwsRegisterUserAlgorithm2(void* sim_data)
{
	/* Register user defined IVP algorithm.�����㷨 */
    MwsIVPSolverProp ivp_prop;
    MwsIVPSolverFcns ivp_fcns;

    ivp_prop.m_name = "myRK45";                        /*���ƣ������ִ�Сд����֤Ψһ*/
    ivp_prop.m_desc = "MYRK45";                        /*�㷨����*/
    ivp_prop.m_fixedStep = moFalse;                      /*�Ƿ�Ϊ��������moTrue��moFalse*/
    ivp_prop.m_ivpType = MWS_IVP_ODE;                   /*�����㷨���ͣ�ODEΪMWS_IVP_ODE��DAEΪMWS_IVP_DAE*/
    ivp_fcns.m_createPtr = &myRK45Create;              /*�����㷨��������ָ��*/
    ivp_fcns.m_createPBPtr = &myRK45ProblemCreate;     /*���������������ָ��*/
    ivp_fcns.m_destroyPBPtr = &myRK45ProblemDestroy;   /*���������������ָ��*/
    ivp_fcns.m_destroyPtr = &myRK45Destroy;            /*�����㷨��������ָ��*/
    ivp_fcns.m_initPtr = &myRK45Init;                  /*��ʼ������ָ��*/
    ivp_fcns.m_interpolatePtr = &myRK45Interpolate;    /*��ֵ��������*/
    ivp_fcns.m_solvePtr = &myRK45Solve;                /*��⺯��ָ��*/
    isimRegisterIVPSolver(sim_data, &ivp_prop, &ivp_fcns);
}

void MwsUnregisterUserAlgorithm1(void* mdl_data)
{
	/* Unregister user defined LS and NLS algorithm. */
}

void MwsUnregisterUserAlgorithm2(void* sim_data)
{
	/* Unregister user defined IVP algorithm. */
    isimUnregisterIVPSolver(sim_data, "myRK45");
}