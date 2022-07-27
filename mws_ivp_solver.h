/***************************************************************************
///
/// Copyright (c) 2020, 苏州同元软控信息技术有限公司
/// All rights reserved.
///
/// @file           mws_ivp_solver.h
/// @brief          积分算法
///
/// @version        v1.0
/// @author         田显钊
/// @date           2020/05/25
///
***************************************************************************/

#ifndef MWS_IVP_SOLVER_H
#define MWS_IVP_SOLVER_H

#include "mws_types.h"

MW_EXTERN_C_BEGIN

/* DAE的形式为:
 *  F(t,y,y')=0
 *  y(t0)=y0
 *  y'(t0)=y'0
 * ODE的形式为:
 *  y'=f(t,y)
 *  y(t0)=y0
 */

typedef void* MwsIVPSolverObj;  /* 积分算法对象 */
typedef void* MwsIVPObj;        /* 求解问题对象 */

/**
 * @brief 积分算法函数返回状态
 */
typedef enum
{
    MWS_IVP_SUCCESS = 0,    /* 成功 */
    MWS_IVP_WARNING = 1,    /* 有警告 */
    MWS_IVP_TSTOP_RETURN,   /* 到达终止时间而返回 */
    MWS_IVP_MEM_FAIL,       /* 内存分配失败 */
    MWS_IVP_FAIL,           /* 失败（原因未知） */
    MWS_IVP_INVALID_INPUT,  /* 输入无效 */
    MWS_IVP_RHSFN_FAIL,     /* ODE通用形式中右端函数出错 */
    MWS_IVP_RESFN_FAIL,     /* DAE残余函数出错 */
}MwsIVPStatus;

/**
 * @brief 算法求解问题类型
 */
typedef enum
{
    MWS_IVP_ODE,    /* ODE问题 */
    MWS_IVP_DAE,    /* DAE问题 */
}MwsIVPType;

/**
 * @brief 工具函数(平台提供，积分算法使用)
 */
typedef struct
{
    /*
     * @breief 日志函数
     * @param[in] user_data     用户数据，创建对象接口的最后一个参数，下同
     * @param[in] error_code    错误号
     * @param[in] where         错误来自何处
     * @param[in] msg           错误信息
     */
    void (*m_logger)(void* user_data, MwsInteger error_code, MwsString where, MwsString msg);

    /*
     * @breief 内存分配函数
     * @param[in] user_data     用户数据
     * @param[in] nobj          对象数量
     * @param[in] size          对象大小（字节数）
     * @return  成功则返回分配的首地址，否则返回空
     */
    void* (*m_allocMemory)(void* user_data, MwsSize nobj, MwsSize size);

    /*
     * @breief 内存释放函数
     * @param[in] user_data     用户数据
     * @param[in] p             内存首地址
     */
    void (*m_freeMemory)(void* user_data, void* p);

    /*
     * @breief 数据内存分配函数
     * @param[in] user_data     用户数据
     * @param[in] nobj          对象数量
     * @param[in] size          对象大小（字节数）
     * @return  成功则返回分配的首地址，否则返回空
     */
    void* (*m_allocDataMemory)(void* user_data, MwsSize nobj, MwsSize size);

    /*
     * @breief 数据内存释放函数
     * @param[in] user_data     用户数据
     * @param[in] p             内存首地址
     */
    void (*m_freeDataMemory)(void* user_data, void* p);

}MwsIVPUtilFcns;

/**
 * @brief 积分算法属性
 */
typedef struct MwsIVPSolverProp
{
    MwsString m_name;           /* 名称（不区分大小写），保证唯一 */
    MwsString m_desc;           /* 算法描述 */
    MwsBoolean m_fixedStep;     /* 是否定步长 */
    MwsIVPType m_ivpType;       /* 算法类型(ODE或DAE) */
}MwsIVPSolverProp;

/**
 * @brief 积分算法选项
 */
typedef struct MwsIVPOptions
{
    MwsBoolean m_stopTimeDefined;       /* 终止时间是否定义 */
    MwsReal m_stopTime;                 /* 终止时间 */
    MwsBoolean m_toleranceDefined;      /* 容许误差是否定义 */
    MwsReal* m_relativeTolerance;       /* 相对容许误差，数组 */
    MwsReal* m_absoluteTolerance;       /* 绝对容许误差，数组 */
    MwsBoolean m_maxStepSizeDefined;    /* 最大步长是否定义 */
    MwsReal m_maxStepSize;              /* 最大步长 */

}MwsIVPOptions;

/*
 * 回调函数（平台提供，积分算法使用）
 */

/*
 * @breief ODE通用形式右端函数 y'=f(t,y)
 * @param[in]user_data 用户数据，创建问题对象接口的最后一个参数，下同
 * @param[in]t          时间
 * @param[in]y          y的值
 * @param[out]yp        y'的值
 * @return 状态，MwsIVPStatus类型的值
 */
typedef MwsInteger (*MwsIVPRshFcnPtr)(void* user_data, MwsReal t, const MwsReal* y, MwsReal* yp);

/*
 * @breief DAE残余量计算回调函数
 * @param[in]user_data 用户数据
 * @param[in]t          时间
 * @param[in]y          y的值
 * @param[in]yp         y'的值
 * @param[out]delta     残余量
 * @return 状态，MwsIVPStatus类型的值
 */
typedef MwsInteger (*MwsIVPResFcnPtr)(void* user_data, MwsReal t, const MwsReal* y, 
    const MwsReal* yp, MwsReal* delta);

/**
 * @brief Jacobian计算回调函数
 * @param[in]user_data  用户数据
 * @param[in]t          时间      
 * @param[in]y          y的值
 * @param[in]yp         y'的值（DAE）
 * @param[in]cj         作用于Jacobian的标量（DAE）
 * @param[out]pd        偏导（pd=dF/dy+cj*dF/dy'或pd=df/dy）    
 * @return 状态，MwsIVPStatus类型的值 
 * @note 
 */
typedef MwsInteger(*MwsIVPJacFcnPtr)(void* user_data, MwsReal t, const MwsReal* y, const MwsReal* yp, 
    MwsReal cj, MwsReal* pd);

/*
 * @breief 积分步完成回调函数
 * @param[in]user_data 用户数据
 * @param[in]t          时间
 * @param[in]y          y的值
 * @return 状态，MwsIVPStatus类型的值
 */
typedef MwsInteger (*MwsIVPStepFinishedPtr)(void* user_data, MwsReal t, const MwsReal* y);

/**
 * @brief 积分算法回调函数
 */
typedef struct  
{
    MwsIVPRshFcnPtr m_rshFunction;          /* ODE右端函数 */
    MwsIVPResFcnPtr m_resFunction;          /* DAE残余函数（与m_rshFunction互斥）*/
    MwsIVPJacFcnPtr m_jacFunction;          /* Jacobian计算函数 */
    MwsIVPStepFinishedPtr m_stepFinished;   /* 积分步完成回调函数 */
}MwsIVPCallback;
 
/*
 * 积分算法接口函数（积分算法提供，平台使用）
 */

/**
  * @breief 积分算法对象创建函数
  * @param[in]util_fcns     工具函数      
  * @param[in]user_data     用户数据(传递给工具函数)
  * @return 积分算法对象
  * @note  工具函数必定会提供
  */
typedef MwsIVPSolverObj (*MwsIVPSolverCreatePtr)(MwsIVPUtilFcns* util_fcns, void* user_data);

/**
 * @brief 求解问题对象创建函数
 * @param[in]solver         积分算法对象
 * @param[in]n              规模，状态变量数量      
 * @param[in]call_back      回调函数
 * @param[in]opt            积分算法选项
 * @param[in]ivp_user_data  用户数据(传递给回调函数)
 * @return 求解问题对象
 * @note  回调函数必定会提供 
 */
typedef MwsIVPObj (*MwsIVPCreatePtr)( MwsIVPSolverObj solver, MwsSize n, MwsIVPCallback* call_back,
    MwsIVPOptions* opt, void* ivp_user_data);

/*
 * @breief 初始化函数
 * @param[in]solver     积分算法对象
 * @param[in]ivp        求解问题对象
 * @param[in]t0         起始时间
 * @param[in]y0         y的初始值（必定不为空）
 * @param[in]yp0        y'的初始值（DAE）
 * @param[in]is_reinit  是否重新初始化（重启）
 * @param[in]reserve    保留参数，暂不使用
 * @return 状态，取MwsIVPStatus的值
 */
typedef MwsInteger (*MwsIVPInitPtr)(MwsIVPSolverObj solver, MwsIVPObj ivp, MwsReal t0, 
    const MwsReal* y0, const MwsReal* yp0, MwsBoolean is_reinit, void* reserve);

/*
 * @breief 求解函数
 * @param[in]solver     积分算法对象
 * @param[in]ivp        求解问题对象
 * @param[in]step_size  步长（定步长/初始化积分步长）
 * @param[in]t          当前时间
 * @param[in]tout       期望的输出时间   
 * @param[out]tret      求解实际到达的时间   
 * @param[out]yret      y的结果值
 * @param[out]ypret     y'的结果值（DAE）
 * @param[in]reserve    保留参数，暂不使用
 * @return 状态，取MwsIVPStatus的值
 */
typedef MwsInteger (*MwsIVPSolvePtr)(MwsIVPSolverObj solver, MwsIVPObj ivp, MwsReal step_size, MwsReal t, 
    MwsReal tout, MwsReal* tret, MwsReal* yret, MwsReal* ypret, void* reserve);

/*
 * @breief 插值函数
 * @param[in]solver     积分算法对象
 * @param[in]ivp        求解问题对象
 * @param[in]tout       期望的输出时间   
 * @param[out]yret      y的结果值
 * @param[in]reserve    保留参数，暂不使用
 * @return 状态，取MwsIVPStatus的值
 */
typedef MwsInteger (*MwsIVPInterpolatePtr)(MwsIVPSolverObj solver, MwsIVPObj ivp, MwsReal tout, 
    MwsReal* yret, void* reserve);
/*
 * @breief 求解问题对象销毁函数
 * @param[in]solver     积分算法对象
 * @param[in]ivp        求解问题对象
 */
typedef void (*MwsIVPDestroyPtr)(MwsIVPSolverObj solver, MwsIVPObj ivp);

/*
 * @breief 积分算法对象销毁函数
 * @param[in]solver     积分算法对象
 */
typedef void (*MwsIVPSolverDestroyPtr)(MwsIVPSolverObj solver);

/**
 * @brief 积分算法接口函数
 */
typedef struct MwsIVPSolverFcns
{
    MwsIVPSolverCreatePtr m_createPtr;      /* 对象创建的函数指针 */
    MwsIVPCreatePtr       m_createPBPtr;    /* 问题对象创建的函数指针 */
    MwsIVPInitPtr m_initPtr;                /* 初始化的函数指针 */
    MwsIVPSolvePtr m_solvePtr;              /* 求解的函数指针 */
    MwsIVPInterpolatePtr m_interpolatePtr;  /* 插值的函数指针 */
    MwsIVPDestroyPtr m_destroyPBPtr;        /* 问题对象销毁的函数指针 */
    MwsIVPSolverDestroyPtr m_destroyPtr;    /* 对象销毁的函数指针 */
}MwsIVPSolverFcns;


/*
 * @breief 注册积分算法
 * @param[in] sim_data 单步层数据
 * @param[in] prop     求解器属性
 * @param[in] fcns     求解器接口函数
 */
MoBoolean isimRegisterIVPSolver(void* sim_data, const MwsIVPSolverProp* prop, const MwsIVPSolverFcns* fcns);

/*
 * @breief 移除注册的积分算法
 * @param[in] sim_data  单步层数据
 * @param[in] name      积分算法名字
 */
MoBoolean isimUnregisterIVPSolver(void* sim_data, MwsString name);

struct MwsLSSolverProp;
struct MwsLSSolverFcns;
/*
 * @breief 注册线性方程组求解算法
 * @param[in] sim_data 单步层数据
 * @param[in] prop     求解器属性
 * @param[in] fcns     求解器接口函数
 */
MoBoolean isimRegisterLSSolver(void* sim_data, const struct MwsLSSolverProp* prop, const struct MwsLSSolverFcns* fcns);

/*
 * @breief 移除注册的线性方程组求解算法
 * @param[in] sim_data  单步层数据
 * @param[in] name      积分算法名字
 */
MoBoolean isimUnregisterLSSolver(void* sim_data, MwsString name);

struct MwsNLSSolverProp;
struct MwsNLSSolverFcns;
/*
 * @breief 注册非线性方程组求解算法
 * @param[in] sim_data 单步层数据
 * @param[in] prop     求解器属性
 * @param[in] fcns     求解器接口函数
 */
MoBoolean isimRegisterNLSSolver(void* sim_data, const struct MwsNLSSolverProp* prop, const struct MwsNLSSolverFcns* fcns);

/*
 * @breief 移除注册的非线性方程组求解算法
 * @param[in] sim_data  单步层数据
 * @param[in] name      积分算法名字
 */
MoBoolean isimUnregisterNLSSolver(void* sim_data, MwsString name);

MW_EXTERN_C_END

#endif /* !MWS_IVP_SOLVER_H */

/***************************************************************************
//   end of file
***************************************************************************/

