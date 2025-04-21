import numpy as np
import tensorflow as tf
import matplotlib.pyplot as plt
from sklearn.feature_selection import mutual_info_regression
import scipy.io
from datetime import datetime
# data
var_h = 0.05

def gen_complex_h():
    real_part = np.random.normal(0., np.sqrt(var_h), [data_size, 1])
    imag_part = np.random.normal(0., np.sqrt(var_h), [data_size, 1])
    return real_part + 1j * imag_part

def gen_complex_n(var_n):
    real_part = np.random.normal(0., np.sqrt(var_n), [data_size, 1])
    imag_part = np.random.normal(0., np.sqrt(var_n), [data_size, 1])
    return real_part + 1j * imag_part

data_size = 300000

snr_arr = np.linspace(0, 30, 31)
H = 128
n_epochs = 300
exper_num = 5
snr_num = len(snr_arr)

class MINE(tf.keras.Model):
    def __init__(self, H):
        super(MINE, self).__init__()
        self.dense1 = tf.keras.layers.Dense(H)
        self.dense2 = tf.keras.layers.Dense(H)
        self.dense3 = tf.keras.layers.Dense(H)
        self.dense4 = tf.keras.layers.Dense(H)
        self.relu = tf.keras.layers.ReLU()
        self.output_layer = tf.keras.layers.Dense(1)

    def call(self, x_conc, y_conc):
        layerx = self.dense1(x_conc)
        layery = self.dense2(y_conc)
        layerx1 = self.relu(layerx)
        layery1 = self.relu(layery)
        layerx2 = self.dense3(layerx1)
        layery2 = self.dense4(layery1)
        layer2 = self.relu(layerx2 + layery2)
        output = self.output_layer(layer2)

        # layerx = self.dense1(x_conc)
        # layery = self.dense2(y_conc)
        # layer2 = self.relu(layerx + layery)
        # output = self.output_layer(layer2)
        return output
        # layerx = self.dense1(x_conc)
        # layery = self.dense2(y_conc)
        # layerx1 = self.relu(layerx)
        # layery1 = self.relu(layery)
        # layerx2 = self.dense3(layerx1)
        # layery2 = self.dense4(layery1)
        # layerx3 = self.relu(layerx2)
        # layery3 = self.relu(layery2)
        # output = self.output_layer(layerx3+layery3)
        # return output

def compute_loss(model, x_in, y_in):
    # shuffle and concatenate
    y_shuffle = tf.random.shuffle(y_in)
    x_conc = tf.concat([x_in, x_in], axis=0)
    y_conc = tf.concat([y_in, y_shuffle], axis=0)

    # propagate the forward pass
    output = model(x_conc, y_conc)

    # split in T_xy and T_x_y predictions
    N_samples = tf.shape(x_in)[0]
    T_xy = output[:N_samples]
    T_x_y = output[N_samples:]

    # compute the negative loss (maximise loss == minimise -loss)
    neg_loss = -(tf.reduce_mean(T_xy) - tf.math.log(tf.reduce_mean(tf.math.exp(T_x_y))))
    return neg_loss

# Create the model and optimizer

# Train the model
MIs = []

MIs_arr = np.zeros((snr_num, exper_num))

mi_numeri = np.zeros(snr_num)

tic = datetime.now()

for snr_idx in range(snr_num):

    snr = snr_arr[snr_idx]
    var_n = (2 * var_h / (10 ** (snr_arr[snr_idx] / 10)))/2
    # ha = 4 * ((2 * var_h) ** 4)
    # hb = ha
    # na = 2 * ((2 * var_h) ** 3) * 2 * var_n + 2 * ((2 * var_h) ** 2) * 2 * var_n + (2 * var_h) * 2 * var_n + 2 * var_n
    # nb = na
    ha = 2*var_h
    hb = ha
    na = 2*var_n
    nb = na

    mi_numeri[snr_idx] = np.log((ha + na) * (hb + nb) / ((ha + na) * (hb + nb) - ha * hb))

    print("SNR = ", snr)
    toc = datetime.now()
    print('Elapsed time: %f seconds' % (toc-tic).total_seconds())

    for exper_idx in range(exper_num):
        model = MINE(H)
        optimizer = tf.keras.optimizers.Adam(learning_rate=0.01)
        for epoch in range(n_epochs):
            # hab1 = gen_complex_h()
            # hba1 = hab1
            # hab2 = gen_complex_h()
            # hba2 = hab2
            # nb11 = gen_complex_n(var_n)
            # na12 = gen_complex_n(var_n)
            # nb23 = gen_complex_n(var_n)
            # na24 = gen_complex_n(var_n)
            # na21 = gen_complex_n(var_n)
            # nb22 = gen_complex_n(var_n)
            # na13 = gen_complex_n(var_n)
            # nb14 = gen_complex_n(var_n)
            # ha = hab1 * hba1 * hab2 * hba2
            # hb = ha
            # na = hba1 * hab2 * hba2 * nb11 + hab2 * hba2 * na12 + hba2 * nb23 + na24
            # nb = hab1 * hba1 * hab2 * na21 + hab1 * hba1 * nb22 + hab1 * na13 + nb14
            ha = gen_complex_h()
            hb = ha
            na = gen_complex_n(var_n)
            nb = gen_complex_n(var_n)
            ha_hat = ha + na
            hb_hat = hb + nb
            ha_hat_stack = np.hstack((ha_hat.real, ha_hat.imag))
            hb_hat_stack = np.hstack((hb_hat.real, hb_hat.imag))

            with tf.GradientTape() as tape:
                # Compute the loss
                neg_loss = compute_loss(model, ha_hat_stack, hb_hat_stack)

            # Compute gradients and update weights
            grads = tape.gradient(neg_loss, model.trainable_variables)
            optimizer.apply_gradients(zip(grads, model.trainable_variables))

            # Save the loss
            MIs.append(-neg_loss.numpy())

        # fig, ax = plt.subplots()
        # ax.plot(range(len(MIs)), MIs, label='MINE estimate')
        # ax.set_xlabel('training steps')
        # ax.legend(loc='best')
        # plt.show()
        # print("MINE estimate", max(MIs))

        MIs_arr[snr_idx][exper_idx] = MIs[-1]


# Plotting
MIs_arr_mean = np.mean(MIs_arr,axis=1)

print("MINE estimate", MIs_arr_mean)
print("Gaussian estimate", mi_numeri)
Gaussian_estimate = mi_numeri

data1 = {'Gaussian_estimate_simple': Gaussian_estimate}
scipy.io.savemat('Gaussian_estimate_simple.mat', data1)

data2 = {'MIs_arr_simple': MIs_arr}
scipy.io.savemat('MIs_arr_simple.mat', data2)

fig, ax = plt.subplots()
ax.plot(snr_arr, MIs_arr_mean, label='MINE estimate')
ax.plot(snr_arr, mi_numeri, label='Theoretical')
# ax.plot([0, len(MIs)], [mi_numeri, mi_numeri], label='True Mutual Information')
# ax.plot([0, len(MIs)], [mi_mutual_info_regression, mi_mutual_info_regression], label='mi_mutual_info_regression')
ax.set_xlabel('SNR')
ax.legend(loc='best')
fig.savefig('MINE_simple.png')
fig.show()