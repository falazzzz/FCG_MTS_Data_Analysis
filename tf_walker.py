# 2018/4/4
# TRY:采用Tensorflow做三元线性回归

import tensorflow as tf
import numpy as np
import read_data
import mts_analysis
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

show = 1
save = 0
sequence = np.array(["yang-baoban_Lu-420-01","yang-baoban_Lu-420-02","yang-baoban_Lu-420-03","yang-baoban_Lu-420-04"])

dadn1, cycles1, dk1 = read_data.ReadMtsResult(sequence=sequence[0])
c1, m1 = mts_analysis.ParisFitting(dadn=dadn1, dk=dk1)
print("Specimen:", sequence[0], ",Paris:c=", str(c1), ",m=", str(m1))
# 读#01数据并拟合Paris

dadn2, cycles2, dk2 = read_data.ReadMtsResult(sequence=sequence[1])
c2, m2 = mts_analysis.ParisFitting(dadn=dadn2, dk=dk2)
print("Specimen:", sequence[1], ",Paris:c=", str(c2), ",m=", str(m2))
# 读#02数据并拟合Paris

threshold = 15
dadn3_temp, cycles3_temp, dk3_temp = read_data.ReadMtsResult(sequence=sequence[2])
dadn3 = dadn3_temp
cycles3 =cycles3_temp
dk3 = dk3_temp

dadn3 = mts_analysis.DataSelectByThreshold(threshold=threshold, parameter=dk3_temp, data=dadn3_temp)
cycles3 = mts_analysis.DataSelectByThreshold(threshold=threshold, parameter=dk3_temp, data=cycles3_temp)
dk3 = mts_analysis.DataSelectByThreshold(threshold=threshold, parameter=dk3_temp, data=dk3_temp)

c3, m3 = mts_analysis.ParisFitting(dadn=dadn3, dk=dk3)
print("Specimen:", sequence[2], ",Paris:c=", str(c3), ",m=", str(m3))
# 读#03数据，筛去DK<15的部分，拟合Paris

dadn4, cycles4, dk4 = read_data.ReadMtsResult(sequence=sequence[3])
c4, m4 = mts_analysis.ParisFitting(dadn=dadn4, dk=dk4)
print("Specimen:", sequence[3], ",Paris:c=", str(c4), ",m=", str(m4))
# 读#04数据并拟合Paris

dadn_r01 = np.concatenate((dadn1, dadn2))
dk_r01 = np.concatenate((dk1, dk2))
c_r01, m_r01 = mts_analysis.ParisFitting(dadn=dadn_r01, dk=dk_r01)
print("r = 0.1:", ",Paris:c=", str(c_r01), ",m=", str(m_r01))
# 合并#01和#02数据(r=0.1)并拟合Paris


dadn_r07 = np.concatenate((dadn3, dadn4))
dk_r07 = np.concatenate((dk3, dk4))
c_r07, m_r07 = mts_analysis.ParisFitting(dadn=dadn_r07, dk=dk_r07)
print("r = 0.7:", ",Paris:c=", str(c_r07), ",m=", str(m_r07))
# 合并#03和#04数据(r=0.7)并拟合Paris


# Tensorflow部分
sess = tf.Session()
# 加载数据集
dk = np.concatenate((dk_r01, dk_r07))
r = np.concatenate((np.full(len(dk_r01), 0.1), np.full(len(dk_r07), 0.7)))
dadn = np.concatenate((dadn_r01, dadn_r07))
logdk = np.log(dk)
logr = np.log(1 - r)
logdadn = np.log(dadn)
x_vals = np.array([[logdk[seq], logr[seq]] for seq in range(len(logdk))])
y_vals = np.array([logdadn[seq] for seq in range(len(logdadn))])

# 申明参数和输出
batch_size = 100
learning_rate = 0.0005
x_data = tf.placeholder(shape=[None, 2], dtype=tf.float32)
y_target = tf.placeholder(shape=[None, 1], dtype=tf.float32)
A = tf.Variable(tf.random_normal(shape=[2, 1]))
b = tf.Variable(tf.random_normal(shape=[1, 1]))
model_output = tf.add(tf.matmul(x_data, A), b)

# 损失函数定义
elastic_param1 = tf.constant(1.)
elastic_param2 = tf.constant(1.)
l1_a_loss = tf.reduce_mean(tf.abs(A))
l2_a_loss = tf.reduce_mean(tf.square(A))
e1_term = tf.multiply(elastic_param1, l1_a_loss)
e2_term = tf.multiply(elastic_param2, l2_a_loss)
loss = tf.expand_dims(tf.add(tf.add(tf.reduce_mean(tf.square(y_target - model_output)), e1_term), e2_term), 0)

# 初始化和训练
init = tf.global_variables_initializer()
sess.run(init)
my_opt = tf.train.GradientDescentOptimizer(learning_rate)
train_step = my_opt.minimize(loss)
loss_vec = []
for i in range(75000):
    rand_index = np.random.choice(len(x_vals), size=batch_size)
    rand_x = x_vals[rand_index]
    rand_y = np.transpose([y_vals[rand_index]])
    sess.run(train_step, feed_dict={x_data: rand_x, y_target: rand_y})
    temp_loss = sess.run(loss, feed_dict={x_data: rand_x, y_target:rand_y})
    loss_vec.append(temp_loss[0])
    if (i+1) % 1000 == 0:
        print('Step:' + str(i+1) + ',A = ' + str(sess.run(A)) + ',b=' + str(sess.run(b)))
        print('Loss=' + str(temp_loss))

# 绘制Loss图像
plt.plot(loss_vec, 'k-')
plt.title('Loss per Generation')
plt.xlabel('Generation')
plt.ylabel('Loss')
plt.show()

c0 = np.exp(sess.run(b))
parameter = sess.run(A)
m0 = parameter[0]
gamma = parameter[1]/m0 + 1
print('WalkerModel Parameter:c0='+str(c0)+',m0='+str(m0)+',gamma='+str(gamma))

# Walker模型计算
dadn_walker = []
dk_walker, r_walker = np.mgrid[10:41:1, 0:0.8:0.1]
dadn_walker = mts_analysis.WalkerCalculating(c0=c0, m0=m0, gamma=gamma, dk=dk_walker, r=r_walker)

# 绘制Walker模型图像
# Plotting: (1)da/dN - dk - r (MTS) plot
r01 = np.full(len(dk_r01), 0.1)
r07 = np.full(len(dk_r07), 0.7)
fig = plt.figure(figsize=(7, 5), dpi=320, num=1)
ax = plt.subplot(111, projection='3d')
ax.scatter(np.log(dk_r01), np.log(1-r01), np.log(dadn_r01), s=1, label='$r=0.1$', color='red', lw=1)
ax.scatter(np.log(dk_r07), np.log(1-r07), np.log(dadn_r07), s=1, label='$r=0.7$', color='blue', lw=1)
Axes3D.plot_surface(self=ax, X=np.log(dk_walker), Y=np.log(1-r_walker), Z=np.log(dadn_walker), rstride=1, cstride=1, cmap=plt.get_cmap('rainbow'))
ax.set_xlabel("DeltaK Applied (MPa*m^0.5)")
ax.set_ylabel("r")
ax.set_zlabel("da/dN (mm/cycle)")
plt.show()
