import pandas as pd
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import roc_curve, auc, log_loss
from sklearn.model_selection import KFold
import joblib
from sklearn.metrics import precision_recall_curve
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
import numpy as np
import os
import tensorflow as tf
import random

# 设置随机种子
tf.random.set_seed(66)
os.environ['PYTHONHASHSEED'] = str(66)
os.environ['TF_DETERMINISTIC_OPS'] = '1'
np.random.seed(66)
random.seed(66)

# 加载数据
data = pd.read_csv('E:/github/new.data/tcga.csv')
# 假设最后一列是目标标签（疾病分类），其他列是基因表达数据
X = data.iloc[:, 1:].values  # 基因表达数据
y = data.iloc[:, 0].values  # 标签（例如，0或1代表两种疾病类型）

# 记录 AUC 和 ROC 曲线
all_fpr = []  # 保存所有的假阳性率
all_tpr = []  # 保存所有的真阳性率
all_roc_auc = []  # 保存所有的 AUC 值

# 新增：记录训练集的 AUC 和 ROC 曲线
all_train_fpr = []
all_train_tpr = []
all_train_roc_auc = []

# 初始化变量来存储当前最佳 AUC
best_auc = 0
best_model = None
best_fold = 0
best_fpr = None
best_tpr = None
best_cutoff = None

# 用于存储每折的训练集和验证集损失及准确率记录
all_train_losses = []
all_val_losses = []
all_train_accuracies = []
all_val_accuracies = []

# 定义十折交叉验证
kf = KFold(n_splits=10, shuffle=True, random_state=42)

for fold, (train_index, test_index) in enumerate(kf.split(X), 1):
    # 划分数据集
    X_train, X_test = X[train_index], X[test_index]
    y_train, y_test = y[train_index], y[test_index]

    # 定义 MLP 模型
    mlp_model = MLPClassifier(hidden_layer_sizes=(100, 50), activation='relu', solver='adam', random_state=42,
                              max_iter=500)

    # 记录训练集和验证集的损失和准确率
    train_losses = []
    val_losses = []
    train_accuracies = []
    val_accuracies = []

    # 自定义训练循环
    for epoch in range(mlp_model.max_iter):
        mlp_model.partial_fit(X_train, y_train, classes=np.unique(y))

        # 计算训练集损失和准确率
        train_loss = mlp_model.loss_
        train_accuracy = mlp_model.score(X_train, y_train)
        train_losses.append(train_loss)
        train_accuracies.append(train_accuracy)

        # 计算验证集损失和准确率
        val_pred_prob = mlp_model.predict_proba(X_test)
        val_loss = log_loss(y_test, val_pred_prob)
        val_accuracy = mlp_model.score(X_test, y_test)
        val_losses.append(val_loss)
        val_accuracies.append(val_accuracy)

    # 测试集预测概率
    y_test_pred_prob = mlp_model.predict_proba(X_test)[:, 1]  # 获取阳性类别的概率

    # 新增：训练集预测概率
    y_train_pred_prob = mlp_model.predict_proba(X_train)[:, 1]

    # 计算概率分布曲线
    kde_0 = gaussian_kde(y_test_pred_prob[y_test == 0])
    kde_1 = gaussian_kde(y_test_pred_prob[y_test == 1])
    x_vals = np.linspace(0, 1, 1000)
    y_vals_0 = kde_0(x_vals)
    y_vals_1 = kde_1(x_vals)

    # 找到相交处的阈值
    diff = np.abs(y_vals_0 - y_vals_1)
    min_index = np.argmin(diff)
    threshold = x_vals[min_index]

    # 获取最佳阈值（基于验证集）
    precision, recall, thresholds = precision_recall_curve(y_test, y_test_pred_prob)
    f1_scores = 2 * (precision * recall) / (precision + recall)
    f1_scores = np.nan_to_num(f1_scores)  # 处理 NaN

    best_threshold = thresholds[np.argmax(f1_scores)]

    # 计算 ROC 曲线
    fpr, tpr, _ = roc_curve(y_test, y_test_pred_prob)
    roc_auc = auc(fpr, tpr)

    # 新增：计算训练集的 ROC 曲线
    train_fpr, train_tpr, _ = roc_curve(y_train, y_train_pred_prob)
    train_roc_auc = auc(train_fpr, train_tpr)

    # 保存当前的假阳性率、真阳性率和 AUC
    all_fpr.append(fpr)
    all_tpr.append(tpr)
    all_roc_auc.append(roc_auc)

    # 新增：保存训练集的假阳性率、真阳性率和 AUC
    all_train_fpr.append(train_fpr)
    all_train_tpr.append(train_tpr)
    all_train_roc_auc.append(train_roc_auc)

    # 保存当前折的训练集和验证集损失及准确率记录
    all_train_losses.append(train_losses)
    all_val_losses.append(val_losses)
    all_train_accuracies.append(train_accuracies)
    all_val_accuracies.append(val_accuracies)

    # 更新最佳模型
    if roc_auc > best_auc:
        best_auc = roc_auc
        best_model = mlp_model
        best_fold = fold
        best_fpr = fpr
        best_tpr = tpr
        best_cutoff = threshold
        # 保存最佳模型
        model_path = f'E:/github/model/best_MLP_model_fold_{best_fold}_auc_{best_auc:.4f}.pkl'
        joblib.dump(best_model, model_path)
        print(f"New best model saved from fold {best_fold} with AUC: {best_auc:.4f}")

    # 保存最佳预测概率和标签用于后续绘图
    np.save('E:/github/model/tissue.best_y_test_pred_prob.npy', y_test_pred_prob)
    np.save('E:/github/model/tissue.best_y_test.npy', y_test)
    np.save('E:/github/model/tissue.best_threshold.npy', best_threshold)

# 打印最佳结果
print(f"\n最佳模型出现在第 {best_fold} 次分组，AUC = {best_auc:.4f}")
print(f"最佳cutoff为: {best_cutoff:.4f}")


# 定义滑动平均平滑函数
def smooth_curve(data, weight=0.7):
    smoothed = []
    last = data[0]
    for point in data:
        smoothed_val = last * weight + (1 - weight) * point
        smoothed.append(smoothed_val)
        last = smoothed_val
    return smoothed


# 对最佳模型所在折的记录数据进行平滑处理
best_train_losses = all_train_losses[best_fold - 1]
best_val_losses = all_val_losses[best_fold - 1]
best_train_accuracies = all_train_accuracies[best_fold - 1]
best_val_accuracies = all_val_accuracies[best_fold - 1]

smooth_train_losses = smooth_curve(best_train_losses)
smooth_val_losses = smooth_curve(best_val_losses)
smooth_train_accuracies = smooth_curve(best_train_accuracies)
smooth_val_accuracies = smooth_curve(best_val_accuracies)

# 绘制平滑后的 Loss 和 Accuracy 曲线
plt.figure(figsize=(14, 6))

# Loss 曲线
plt.subplot(1, 2, 1)
plt.plot(smooth_train_losses, label='Training Loss', color='blue')
plt.plot(smooth_val_losses, label='Validation Loss', color='orange')
plt.xlabel('Epoch')
plt.ylabel('Loss')
plt.title(f'Loss Curve for Best Model (Fold {best_fold})')
plt.legend()

# Accuracy 曲线
plt.subplot(1, 2, 2)
plt.plot(smooth_train_accuracies, label='Training Accuracy', color='blue')
plt.plot(smooth_val_accuracies, label='Validation Accuracy', color='orange')
plt.xlabel('Epoch')
plt.ylabel('Accuracy')
plt.title(f'Accuracy Curve for Best Model (Fold {best_fold})')
plt.legend()

plt.tight_layout()
plt.show()

# 绘制最佳模型的损失曲线
plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
plt.plot(best_train_losses, label='Training Loss')
plt.plot(best_val_losses, label='Validation Loss')
plt.xlabel('Epoch')
plt.ylabel('Loss')
plt.title(f'Loss Curve for Best Model (Fold {best_fold})')
plt.legend()

# 绘制最佳模型的准确率曲线
plt.subplot(1, 2, 2)
plt.plot(best_train_accuracies, label='Training Accuracy')
plt.plot(best_val_accuracies, label='Validation Accuracy')
plt.xlabel('Epoch')
plt.ylabel('Accuracy')
plt.title(f'Accuracy Curve for Best Model (Fold {best_fold})')
plt.legend()

plt.tight_layout()
plt.show()

# 绘制最佳模型的 ROC 曲线
plt.figure(figsize=(8, 6))
plt.plot(best_fpr, best_tpr, label=f'Best ROC Curve (AUC = {best_auc:.4f})', color='blue', lw=2)
plt.plot([0, 1], [0, 1], 'k--', label='Random Guess', lw=1)
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title(f'ROC Curve for Best Model (Fold {best_fold})')
plt.legend(loc='lower right')
plt.show()

# 绘制所有折的训练集的 ROC 曲线
plt.figure(figsize=(10, 8))
for i in range(len(all_fpr)):
    plt.plot(all_fpr[i], all_tpr[i], lw=2, alpha=0.7, label=f'Train ROC fold {i + 1} (AUC = {all_roc_auc[i]:.2f})')

plt.plot([0, 1], [0, 1], color='navy', linestyle='--', lw=2)
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Train Dataset ROC Curves')
plt.legend(loc="lower right", fontsize='small', ncol=2)
plt.show()

# 新增：绘制最佳模型下训练集的 ROC 曲线
best_train_fpr = all_train_fpr[best_fold - 1]
best_train_tpr = all_train_tpr[best_fold - 1]
best_train_roc_auc = all_train_roc_auc[best_fold - 1]

plt.figure(figsize=(8, 6))
plt.plot(best_train_fpr, best_train_tpr, label=f'Best Training ROC Curve (AUC = {best_train_roc_auc:.4f})',
         color='green', lw=2)
plt.plot([0, 1], [0, 1], 'k--', label='Random Guess', lw=1)
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title(f'Training ROC Curve for Best Model (Fold {best_fold})')
plt.legend(loc='lower right')
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
