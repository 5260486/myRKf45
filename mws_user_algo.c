#include "mws_ivp_solver.h"
#include "mws_ls_solver.h"
#include "mws_nls_solver.h"

#include "my_euler.c"    /*�Զ����㷨ͷ�ļ�*/

void MwsRegisterUserAlgorithm1(void* mdl_data)
{
	/* Register user defined LS and NLS algorithm. ���Ի�������㷨*/
}

void MwsRegisterUserAlgorithm2(void* sim_data)
{
	/* Register user defined IVP algorithm.�����㷨 */
    MwsIVPSolverProp ivp_prop;
    MwsIVPSolverFcns ivp_fcns;

    ivp_prop.m_name = "myeuler";                        /*���ƣ������ִ�Сд����֤Ψһ*/
    ivp_prop.m_desc = "MYEULER";                        /*�㷨����*/
    ivp_prop.m_fixedStep = moTrue;                      /*�Ƿ�Ϊ��������moTrue��moFalse*/
    ivp_prop.m_ivpType = MWS_IVP_ODE;                   /*�����㷨���ͣ�ODEΪMWS_IVP_ODE��DAEΪMWS_IVP_DAE*/
    ivp_fcns.m_createPtr = &myEulerCreate;              /*�����㷨��������ָ��*/
    ivp_fcns.m_createPBPtr = &myEulerProblemCreate;     /*���������������ָ��*/
    ivp_fcns.m_destroyPBPtr = &myEulerProblemDestroy;   /*���������������ָ��*/
    ivp_fcns.m_destroyPtr = &myEulerDestroy;            /*�����㷨��������ָ��*/
    ivp_fcns.m_initPtr = &myEulerInit;                  /*��ʼ������ָ��*/
    ivp_fcns.m_interpolatePtr = &myEulerInterpolate;    /*��ֵ��������*/
    ivp_fcns.m_solvePtr = &myEulerSolve;                /*��⺯��ָ��*/
    isimRegisterIVPSolver(sim_data, &ivp_prop, &ivp_fcns);
}

void MwsUnregisterUserAlgorithm1(void* mdl_data)
{
	/* Unregister user defined LS and NLS algorithm. */
}

void MwsUnregisterUserAlgorithm2(void* sim_data)
{
	/* Unregister user defined IVP algorithm. */
    isimUnregisterIVPSolver(sim_data, "myeuler");
}