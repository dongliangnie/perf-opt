项目：卡尔曼滤波(KF)算子性能优化

1. 编译环境
   - 平台：鲲鹏 ARMv8 + Linux
   - 编译器：gcc 9.3+ / clang

2. 编译与执行
   编译：gcc -O2 kf_linux.c kf_linux_test.c timestamp.c -o kf_run -lm
   执行：./kf_run

3. 接口说明
   - kf_linux_iopointer: 调用者提供输入输出缓冲区指针
   - kf_linux_ioself_profiling: 自产生数据，计时并打印
   - kf_linux_ioself: 自产生数据，不计时不打印

4. 优化方法
   - 指针步长消减：减少重复索引计算
   - i-l-j循环翻转：优化Cache空间局部性
   - 连续内存访问：B矩阵按行访问，提高Cache命中率
   - 矩阵运算合并：减少临时变量和内存拷贝


5. 文件清单
   - kf_linux.h      : 接口声明
   - kf_linux.c      : 接口实现
   - kf_linux_test.c : 测试主函数
   - timestamp.h     : 时间函数声明
   - timestamp.c     : 时间函数实现
   - readme.txt      : 本文件
