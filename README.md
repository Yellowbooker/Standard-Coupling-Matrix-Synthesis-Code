# Standard-Coupling-Matrix-Synthesis-Code
Standard Coupling Matrix Synthesis Code  
标准化耦合矩阵综合代码：包含广义切比雪夫滤波器耦合矩阵综合与矩阵旋转（折叠型、箭型、级联CTCQ型）代码以及基于广义切比雪夫滤波器的有耗滤波器耦合矩阵综合代码  
请运行主函数：  
GChebyshevFilterMain.m（无耗广义切比雪夫滤波器耦合矩阵综合）  
LossyFilterMain.m（基于广义切比雪夫滤波器的有耗滤波器耦合矩阵综合）  
子函数说明：  
1. General_Chebyshev.m（广义切比雪夫滤波函数的综合）
2. Cheby2EPF.m（广义切比雪夫滤波器对应特征多项式的综合）
3. EPF2Y.m（通过特征多项式E、P、F求解Y参数，无耗情况）
4. EPF2Ylossy.m（通过特征多项式E、P、F求解Y参数，有耗情况）
5. CheckUnitary.m（确认S参数是否满足无耗条件，即幺正性）
6. Y2CMtrans.m（Y参数到横向耦合矩阵）
7. Rotate.m（对耦合矩阵进行一次旋转，消去特定位置的矩阵元素）
8. to_foldedCM.m（旋转至折叠型耦合矩阵，顺次耦合为正，模值小于0.0001的元素置零）
9. CM2arrow.m（旋转至箭型耦合矩阵）
10. TriExtract.m（从箭型耦合矩阵中提取出CT结构，两个CT结构可转换为一个CQ结构）
11. CMFC_Response.m（计算耦合矩阵对应的S参数）
12. LossyCM2PCP.m（有耗滤波器横向耦合矩阵至并联双耦合谐振器耦合矩阵）
13. LossDistribution（有耗滤波器并联双耦合谐振器耦合矩阵的损耗重新分配，双耦合谐振器的Q值一致）
14. residue2.m
