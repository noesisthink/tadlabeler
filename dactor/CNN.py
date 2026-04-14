import os
import numpy as np
from sklearn.metrics import roc_curve, auc, precision_recall_fscore_support, accuracy_score, matthews_corrcoef
from sklearn.model_selection import StratifiedKFold
import tensorflow as tf
from tensorflow.keras.models import Model, load_model
from tensorflow.keras.layers import (
    Dense, Activation, Dropout, Conv1D, MaxPooling1D,
    GlobalMaxPooling1D, Input
)
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.callbacks import EarlyStopping, ModelCheckpoint

# GPU 设置
os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
os.environ["CUDA_VISIBLE_DEVICES"] = "1"

gpus = tf.config.experimental.list_physical_devices('GPU')
if gpus:
    for gpu in gpus:
        tf.config.experimental.set_memory_growth(gpu, True)

# 载入数据
x_train = np.load('GM12878_train_encoded_data.npy')
x_test = np.load('GM12878_test_encoded_data.npy')
x_val = np.load('GM12878_val_encoded_data.npy')
y_train = np.load('GM12878_train_encoded_labels.npy')
y_test = np.load('GM12878_test_encoded_labels.npy')
y_val = np.load('GM12878_val_encoded_labels.npy')

print("DNA序列数据形状:")
print(f"x_train: {x_train.shape}, x_val: {x_val.shape}, x_test: {x_test.shape}")
print("标签形状:")
print(f"y_train: {y_train.shape}, y_val: {y_val.shape}, y_test: {y_test.shape}")

# 模型参数
input_shape = x_train.shape[1:3]
kernel_size = 9
learning_rate = 0.001
num_kernels = 64
name = "GM12878_CNN_manual"
NB_EPOCH = 150
BATCH_SIZE = 16
KERNEL_INITIAL = 'glorot_uniform'
LOSS = 'binary_crossentropy'
METRICS = ['accuracy']

# 保存训练过程
def SaveHistory(Tuning, outfile):
    Hist = np.empty((len(Tuning.history['loss']), 4))
    Hist[:, 0] = Tuning.history['val_loss']
    Hist[:, 1] = Tuning.history['val_accuracy']
    Hist[:, 2] = Tuning.history['loss']
    Hist[:, 3] = Tuning.history['accuracy']
    np.savetxt(outfile, Hist, fmt='%.8f', delimiter=",", header="val_loss,val_accuracy,train_loss,train_acc", comments="")
    return Hist

# 评估
def GetMetrics(model, x, y):
    pred_p = model.predict(x)
    pred = (pred_p > 0.5).astype("int32")
    fpr, tpr, _ = roc_curve(y, pred_p)
    aucv = auc(fpr, tpr)
    precision, recall, fscore, support = precision_recall_fscore_support(y, pred, average='macro', zero_division=0)
    print('AUC:', aucv)
    print('Accuracy:', accuracy_score(y, pred))
    print('MCC:', matthews_corrcoef(y, pred))
    print('Precision:', precision, 'Recall:', recall, 'F1:', fscore)
    return [aucv, accuracy_score(y, pred), matthews_corrcoef(y, pred), precision, recall, fscore, support]

# 构建 CNN 模型（无 LSTM）
def build_CNN_model(INPUT_SHAPE, NUM_KERNEL, KERNEL_SIZE, name):
    dna_input = Input(shape=INPUT_SHAPE, name="DNA_Input")
    x = Conv1D(NUM_KERNEL, kernel_size=KERNEL_SIZE, kernel_initializer=KERNEL_INITIAL)(dna_input)
    x = Activation("relu")(x)
    x = Dropout(0.2)(x)
    x = MaxPooling1D(pool_size=2)(x)

    x = Conv1D(NUM_KERNEL, kernel_size=KERNEL_SIZE, kernel_initializer=KERNEL_INITIAL)(x)
    x = Activation("relu")(x)
    x = Dropout(0.3)(x)
    x = MaxPooling1D(pool_size=2)(x)

    x = Conv1D(NUM_KERNEL, kernel_size=KERNEL_SIZE, kernel_initializer=KERNEL_INITIAL)(x)
    x = Activation("relu")(x)
    x = Dropout(0.3)(x)
    x = MaxPooling1D(pool_size=2)(x)

    x = GlobalMaxPooling1D()(x)
    output = Dense(1, activation="sigmoid")(x)

    model = Model(inputs=dna_input, outputs=output)
    model.compile(loss=LOSS, optimizer=Adam(learning_rate=learning_rate), metrics=METRICS)
    model.summary()
    return model

# 训练和评估
def train_model():
    model = build_CNN_model(input_shape, num_kernels, kernel_size, name)
    filepath = name + "_model.hdf5"
    checkpoint = ModelCheckpoint(filepath, monitor='val_loss', save_best_only=True)
    early_stop = EarlyStopping(monitor='val_loss', patience=8)

    Tuning = model.fit(
        x_train, y_train,
        validation_data=(x_val, y_val),
        batch_size=BATCH_SIZE,
        epochs=NB_EPOCH,
        verbose=1,
        callbacks=[checkpoint, early_stop]
    )

    print("train,", filepath)
    saved_model = load_model(filepath)
    GetMetrics(saved_model, x_train, y_train)

    print("test,", filepath)
    GetMetrics(saved_model, x_test, y_test)

    SaveHistory(Tuning, name + "_training_history.txt")

# 主程序入口
if __name__ == "__main__":
    train_model()
