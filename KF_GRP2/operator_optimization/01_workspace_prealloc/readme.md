阶段 1：workspace 预分配优化

一、阶段定位

本阶段在第 0 阶段“老师要求的 6 文件基础版”基础上进行第一个算子代码优化。
本阶段仍然使用 -O0 编译，目的是观察代码结构优化本身带来的效果，不让 GCC -O2 自动优化干扰判断。

对应算子优化表第 1 行：
    优化前：526.086937 ms
    优化后：524.148291 ms
    加速：0.37%
    final_sum：33407.632487

二、本阶段包含的 6 个文件

1. kf_linux.h
   在第 0 阶段基础上新增 KfWorkspace 结构体。

2. kf_linux.c
   将核心函数内部申请的工作矩阵移动到 KfWorkspace 中，由外层提前分配和初始化。

3. kf_linux_test.c
   测试主函数，仍然测试 horizon=10 到 298。

4. timestamp.h
   时间戳接口声明。

5. timestamp.c
   时间戳接口实现。

6. readme.txt
   本阶段说明。

三、本阶段做了什么优化

第 0 阶段中，kf_linux_iopointer() 内部会申请和释放大量数组：

    x
    P
    Q
    R
    K
    identity
    temp_p_plus_r
    inv_temp_p_plus_r
    temp_i_minus_k
    p_old
    z_minus_x
    k_mult_z_minus_x

这些数组都是 KF 核心计算过程中需要反复使用的工作区。
第 0 阶段的做法是：每次调用核心函数时，在 kf_linux_iopointer() 内部 calloc，算完后 free。

本阶段改成：
    1. 在 kf_linux.h 中定义 KfWorkspace。
    2. 在 KfInput 中加入 KfWorkspace *workspace。
    3. 在 prepare_kf_input() 中调用 create_workspace() 提前分配并初始化工作区。
    4. kf_linux_iopointer() 不再申请和释放这些工作数组，只直接使用 input->workspace。
    5. 结束时由 release_kf_input() 统一释放 workspace。

四、为什么要这么优化

动态内存分配 calloc/free 本身有额外开销。
如果把内存申请放在核心算子函数内部，核心函数既要负责计算，又要负责内存管理，代码结构比较混乱。

workspace 预分配的思想是：
    核心函数只做核心计算；
    输入准备阶段负责内存申请；
    释放阶段统一回收内存。

这样做的好处：
    1. 减少 kf_linux_iopointer() 内部的内存管理开销。
    2. 让核心算子代码更清晰。
    3. 为后续矩阵乘法、循环展开、ARM 对齐等优化打基础。

五、为什么本阶段加速不大

本阶段实测只加速约 0.37%。
这是正常现象。

原因是：
    horizon=298 时，程序主要耗时集中在 mat_mul() 的三层矩阵乘法。
    相比 O(N^3) 的矩阵乘法，calloc/free 不是最大瓶颈。

所以本阶段虽然加速不大，但它是合理的代码结构优化，也为后面优化做铺垫。

六、修改了哪些代码

主要修改 kf_linux.h：

    新增 KfWorkspace：
        double *x;
        double *p;
        double *q;
        double *r;
        double *k;
        double *identity;
        double *temp_p_plus_r;
        double *inv_temp_p_plus_r;
        double *temp_i_minus_k;
        double *p_old;
        double *z_minus_x;
        double *k_mult_z_minus_x;

    KfInput 中新增：
        KfWorkspace *workspace;

主要修改 kf_linux.c：

    新增函数：
        create_workspace()
        init_workspace()
        release_workspace()

    修改 prepare_kf_input()：
        生成 measurements 后，同时创建 workspace。

    修改 kf_linux_iopointer()：
        删除内部 calloc/free 工作区的逻辑，改为直接使用 workspace。

七、本阶段没有改变什么

本阶段没有改变：
    KF 数学公式
    dense 完整矩阵计算方式
    mat_mul() 的 i-j-l 原始循环顺序
    mat_add() 和 mat_sub() 的普通循环
    P 矩阵更新方式
    编译选项 -O0
    ARM 平台专项优化

八、优化类型

本阶段属于：
    KF 算子代码优化
    内存管理优化
    工作区预分配优化
    工程结构优化

九、编译方法

    gcc -Wall -Wextra -O0 kf_linux.c kf_linux_test.c timestamp.c -lm -o kf_step1_workspace_o0

十、运行方法

    ./kf_step1_workspace_o0 | tee step1_workspace_o0_result.txt

十一、正确性检查

horizon=298 时，参考输出：
    time_ms 约为 524.148291
    final_sum 应为 33407.632487

如果 final_sum 与第 0 阶段一致，说明 workspace 预分配没有改变 KF 计算结果。
