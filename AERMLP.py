import pandas as pd
import numpy as np
from sklearn.model_selection import KFold
from tensorflow.keras.models import Model, Sequential
from tensorflow.keras.layers import Input, Dense, GRU, BatchNormalization, Dropout
from sklearn.metrics import  precision_recall_curve
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.regularizers import l2
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt
from tqdm import tqdm
from tensorflow.keras.callbacks import Callback
from imblearn.over_sampling import SMOTE
import pickle
import os
import tensorflow as tf
import random

# 设置随机种子
tf.random.set_seed(66)
os.environ['PYTHONHASHSEED'] = str(66)
os.environ['TF_DETERMINISTIC_OPS'] = '1'
np.random.seed(66)
random.seed(66)

# 数据加载
data = pd.read_csv('E:/github/new.data/tcga.csv')
X = data.iloc[:, 1:].values  # 特征数据
y = data.iloc[:, 0].values  # 标签数据

# 保存结果
all_fpr = []
all_tpr = []
all_roc_auc = []
all_train_fpr = []
all_train_tpr = []
all_train_roc_auc = []

best_auc = 0
best_model = None
best_fold = 0

best_thresholds = []

# 定义 10 折交叉验证
kf = KFold(n_splits=10, shuffle=True, random_state=42)

# 保存每一折的数据
split_data = []

# 开始 10 折交叉验证
for fold, (train_index, test_index) in enumerate(kf.split(X), 1):
    # 划分数据集
    X_train, X_test = X[train_index], X[test_index]
    y_train, y_test = y[train_index], y[test_index]

    # 使用 SMOTE 进行过采样
    smote = SMOTE(random_state=42)
    X_train, y_train = smote.fit_resample(X_train, y_train)

    batch_count = int(np.ceil(len(X_train) / 64))
    print(f"每个 epoch 的 batch 数量: {batch_count}")

    split_data.append({
        "X_train": X_train.tolist(),
        "X_test": X_test.tolist(),
        "y_train": y_train.tolist(),
        "y_test": y_test.tolist()})

    # 自动编码器部分
    input_dim = X_train.shape[1]
    encoding_dim = 64  # 编码器输出的维度

    input_layer = Input(shape=(input_dim,))
    encoded = Dense(128, activation='relu')(input_layer)
    encoded = Dense(encoding_dim, activation='relu')(encoded)
    decoded = Dense(128, activation='relu')(encoded)
    decoded = Dense(input_dim, activation='linear')(decoded)

    autoencoder = Model(input_layer, decoded)
    encoder = Model(input_layer, encoded)

    autoencoder.compile(optimizer=Adam(learning_rate=0.0005), loss='mse')
    autoencoder.fit(X_train, X_train, epochs=100, batch_size=64, verbose=0)

    # 使用编码器对数据进行编码
    X_train_encoded = encoder.predict(X_train)
    X_test_encoded = encoder.predict(X_test)

    # 调整数据维度为 (样本数, 时间步长, 特征数)
    X_train_encoded = X_train_encoded.reshape(X_train_encoded.shape[0], 1, X_train_encoded.shape[1])  # 时间步长为 1
    X_test_encoded = X_test_encoded.reshape(X_test_encoded.shape[0], 1, X_test_encoded.shape[1])

    sequence_length = X_train_encoded.shape[1]
    input_dim_encoded = X_train_encoded.shape[2]

    # 定义模型
    model = Sequential()
    model.add(GRU(128, input_shape=(sequence_length, input_dim_encoded), return_sequences=True, activation='tanh', kernel_regularizer=l2(0.01)))
    model.add(BatchNormalization())
    model.add(Dropout(0.5))
    model.add(GRU(64, return_sequences=False, activation='tanh', kernel_regularizer=l2(0.01)))
    model.add(BatchNormalization())
    model.add(Dropout(0.5))

    model.add(Dense(64, activation='gelu'))  # MLP 第一层
    model.add(BatchNormalization())
    model.add(Dropout(0.5))
    model.add(Dense(32, activation='gelu'))  # MLP 第二层
    model.add(BatchNormalization())
    model.add(Dropout(0.5))
    model.add(Dense(1, activation='sigmoid'))  # 输出层

    model.compile(optimizer=Adam(learning_rate=0.0001), loss='binary_crossentropy', metrics=['accuracy'])

    # 定义进度条
    class TQDMProgressBar(Callback):
        def on_epoch_begin(self, epoch, logs=None):
            if epoch == 0:
                self.epochs_bar = tqdm(total=self.params['epochs'], desc="Training Progress", ncols=100)

        def on_epoch_end(self, epoch, logs=None):
            self.epochs_bar.update(1)
            self.epochs_bar.set_postfix(loss=logs['loss'], accuracy=logs['accuracy'])
            if epoch == self.params['epochs'] - 1:
                self.epochs_bar.close()

    progress_bar = TQDMProgressBar()

    # 模型训练
    history = model.fit(
        X_train_encoded,  # 使用编码后的数据
        y_train,
        epochs=3000,
        batch_size=64,
        verbose=0,
        validation_data=(X_test_encoded, y_test),
        callbacks=[progress_bar]
    )

    # 测试集预测
    y_test_pred_prob = model.predict(X_test_encoded)
    fpr, tpr, thresholds = roc_curve(y_test, y_test_pred_prob)
    roc_auc = auc(fpr, tpr)

    # 获取最佳阈值（基于验证集）
    precision, recall, thresholds = precision_recall_curve(y_test, y_test_pred_prob)
    f1_scores = 2 * (precision * recall) / (precision + recall)
    f1_scores = np.nan_to_num(f1_scores)  # 处理 NaN

    best_threshold = thresholds[np.argmax(f1_scores)]

    # 训练集预测
    y_train_pred_prob = model.predict(X_train_encoded)
    train_fpr, train_tpr, _ = roc_curve(y_train, y_train_pred_prob)
    train_roc_auc = auc(train_fpr, train_tpr)

    # 保存结果
    all_fpr.append(fpr)
    all_tpr.append(tpr)
    all_roc_auc.append(roc_auc)
    all_train_fpr.append(train_fpr)
    all_train_tpr.append(train_tpr)
    all_train_roc_auc.append(train_roc_auc)
    print(f"Fold {fold} Test AUC: {roc_auc:.4f}, Train AUC: {train_roc_auc:.4f}")

    # 如果当前 AUC 值比最佳 AUC 更高，则保存模型和自动编码器权重
    if roc_auc > best_auc:
        best_auc = roc_auc
        best_fold = fold
        # 保存当前最好的模型
        model.save(f'E:/github/model/gelu.0.0001AERMLP_best_model_fold_{best_fold}_auc_{best_auc:.4f}.h5')
        # 保存最佳折对应的自动编码器权重
        autoencoder.save_weights(f'E:/github/model/gelu.0.0001autoencoder_weights_best_fold_{best_fold}.weights.h5')
        print(f"New best model saved with AUC: {best_auc:.4f}")

        print(f"New best model saved with AUC: {best_auc:.4f} and threshold: {best_threshold:.4f}")

        import numpy as np

        # 保存最佳预测概率和标签用于后续绘图
        np.save('E:/github/model/best_y_test_pred_prob.npy', y_test_pred_prob)
        np.save('E:/github/model/best_y_test.npy', y_test)
        np.save('E:/github/model/best_threshold.npy', best_threshold)

        # 绘制训练集 ROC 曲线
        plt.figure(figsize=(8, 6))
        plt.plot(train_fpr, train_tpr, label=f'Train ROC Curve (AUC = {train_roc_auc:.4f})', color='green')
        plt.plot([0, 1], [0, 1], 'k--', label='Random Guess')
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title(f'Train ROC Curve for Fold {fold}')
        plt.legend(loc='lower right')

        # 绘制测试集 ROC 曲线
        plt.figure(figsize=(8, 6))
        plt.plot(fpr, tpr, label=f'Test ROC Curve (AUC = {roc_auc:.4f})', color='blue')
        plt.plot([0, 1], [0, 1], 'k--', label='Random Guess')
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title(f'Validation ROC Curve for Fold {fold}')
        plt.legend(loc='lower right')

        # 绘制并保存训练和验证集的 Loss 和 Accuracy 曲线
        plt.figure(figsize=(12, 5))
        # Loss 曲线
        plt.subplot(1, 2, 1)
        plt.plot(history.history['loss'], label='Train Loss')
        plt.plot(history.history['val_loss'], label='Validation Loss')
        plt.xlabel('Epochs')
        plt.ylabel('Loss')
        plt.title('Train and Validation Loss')
        plt.legend()

        # Accuracy 曲线
        plt.subplot(1, 2, 2)
        plt.plot(history.history['accuracy'], label='Train Accuracy')
        plt.plot(history.history['val_accuracy'], label='Validation Accuracy')
        plt.xlabel('Epochs')
        plt.ylabel('Accuracy')
        plt.title('Train and Validation Accuracy')
        plt.legend()

# 保存到二进制文件
#with open("E:/github/model/gelu.0.0001.AERMLP_split_data.pkl", "wb") as f:
    #pickle.dump(split_data, f)

# 绘制所有 ROC 曲线
plt.figure(figsize=(10, 8))
for i in range(len(all_fpr)):
    plt.plot(all_fpr[i], all_tpr[i], lw=2, alpha=0.7, label=f'ROC fold {i + 1} (AUC = {all_roc_auc[i]:.2f})')

# 绘制参考线
plt.plot([0, 1], [0, 1], color='navy', linestyle='--', lw=2)
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Validation Dataset ROC Curves')
plt.legend(loc="lower right", fontsize='small', ncol=2)
plt.show()

# 绘制训练集的 ROC 曲线
plt.figure(figsize=(10, 8))
for i in range(len(all_train_fpr)):
    plt.plot(all_train_fpr[i], all_train_tpr[i], lw=2, alpha=0.7, label=f'Train ROC fold {i + 1} (AUC = {all_train_roc_auc[i]:.2f})')

plt.plot([0, 1], [0, 1], color='navy', linestyle='--', lw=2)
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Train Dataset ROC Curves')
plt.legend(loc="lower right", fontsize='small', ncol=2)
plt.show()

# 绘制 AUC 散点图
import matplotlib.patches as patches
plt.figure(figsize=(8, 6))
# 获取当前的轴对象
ax = plt.gca()
plt.scatter(range(1, len(all_roc_auc) + 1), all_roc_auc, color='blue', marker='o')\
# 添加框，包围所有散点
min_auc = min(all_roc_auc)
max_auc = max(all_roc_auc)
plt.gca().add_patch(patches.Rectangle((0.5, min_auc - 0.02), len(all_roc_auc), max_auc - min_auc + 0.04,
                                     linewidth=1, edgecolor='black', facecolor='none', linestyle='-'))
# 设置纵坐标范围为 0 到 1
plt.ylim(0, 1)
plt.xlim(0, 20)
# 关闭网格线
ax.grid(False)
plt.xlabel('Fold Number')
plt.ylabel('AUC')
plt.title('AUC for Each Fold')
plt.show()

# 假设 history 是训练后的历史对象
losses = history.history['loss']  # 每个 epoch 的训练损失
val_losses = history.history['val_loss']  # 每个 epoch 的验证损失

# 计算每个 epoch 到当前 epoch 的中位数损失
mean_losses = [np.median(losses[:i+1]) for i in range(len(losses))]
mean_val_losses = [np.median(val_losses[:i+1]) for i in range(len(val_losses))]

# 绘制平均损失曲线
epochs = range(1, len(losses) + 1)
plt.plot(epochs, mean_losses, label='Mean Training Loss')
plt.plot(epochs, mean_val_losses, label='Mean Validation Loss')
plt.xlabel('Epochs')
plt.ylabel('Mean Loss')
plt.legend()
plt.show()

# 假设 'accuracy' 和 'val_accuracy' 存储了每个 epoch 的训练和验证准确率
accuracy = history.history['accuracy']
val_accuracy = history.history['val_accuracy']

# 计算每个 epoch 到当前 epoch 的平均准确率
mean_accuracy = [np.mean(accuracy[:i+1]) for i in range(len(accuracy))]
mean_val_accuracy = [np.mean(val_accuracy[:i+1]) for i in range(len(val_accuracy))]

# 绘制平均准确率曲线
epochs = range(1, len(accuracy) + 1)
plt.plot(epochs, mean_accuracy, label='Mean Training Accuracy')
plt.plot(epochs, mean_val_accuracy, label='Mean Validation Accuracy')
plt.xlabel('Epochs')
plt.ylabel('Mean Accuracy')
plt.legend()
plt.show()

# 计算并打印测试集AUC的均值和标准差
all_roc_auc_np = np.array(all_roc_auc)
mean_auc = np.mean(all_roc_auc_np)
std_auc = np.std(all_roc_auc_np)
print(f"测试集AUC的均值: {mean_auc:.4f}")
print(f"测试集AUC的标准差: {std_auc:.4f}")

# 计算并打印训练集AUC的均值和标准差
all_train_roc_auc_np = np.array(all_train_roc_auc)
mean_train_auc = np.mean(all_train_roc_auc_np)
std_train_auc = np.std(all_train_roc_auc_np)
print(f"训练集AUC的均值: {mean_train_auc:.4f}")
print(f"训练集AUC的标准差: {std_train_auc:.4f}")
