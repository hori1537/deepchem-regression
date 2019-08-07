

import os
import tkinter
import tkinter.filedialog
from tkinter import ttk
from tkinter import N,E,S,W
from tkinter import font

import pubchempy as pcp

from mordred import descriptors, Calculator

from rdkit import rdBase, Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import AllChem, Draw, PandasTools

import random

import pandas as pd
import numpy as np


import sklearn
from sklearn.preprocessing import StandardScaler
from sklearn import model_selection
from keras.models import Sequential
from keras.layers import Input, Activation, Conv2D, BatchNormalization, Flatten, Dense
from keras.optimizers import SGD

import sascorer


# refer https://horomary.hatenablog.com/entry/2018/10/21/122025
# refer https://www.ag.kagawa-u.ac.jp/charlesy/2017/07/27/deepchem%E3%81%AB%E3%82%88%E3%82%8B%E6%BA%B6%E8%A7%A3%E5%BA%A6%E4%BA%88%E6%B8%AC-graph-convolution-%E3%83%8B%E3%83%A5%E3%83%BC%E3%83%A9%E3%83%AB%E3%83%8D%E3%83%83%E3%83%88%E3%83%AF%E3%83%BC%E3%82%AF/
# refer https://future-chem.com/rdkit-intro/ 
#モジュールの読み込み
import tensorflow as tf
import deepchem as dc
import numpy as np
from deepchem.models.tensorgraph.models.graph_models import GraphConvModel
import sys
from sys import exit



def get_csv():
    tk_c = tkinter.Tk()
    current_dir = os.getcwd()

    csv_file_path = tkinter.filedialog.askopenfilename(initialdir = current_dir, 
    title = 'choose the csv', filetypes = [('csv file', '*.csv')])

    t_csv.set(csv_file_path)

    


def graph_conv_training():
    graph_featurizer = dc.feat.graph_features.ConvMolFeaturizer()

    loader = dc.data.data_loader.CSVLoader( tasks = [t_task.get()], smiles_field = t_smiles.get(), id_field = t_id.get(), featurizer = graph_featurizer )
    dataset = loader.featurize( t_csv.get() )

    splitter = dc.splits.splitters.RandomSplitter()
    trainset, testset = splitter.train_test_split( dataset )
    
    hp = dc.molnet.preset_hyper_parameters
    param = hp.hps[ 'graphconvreg' ]
    print( param )
    
    batch_size = 48

    from deepchem.models.tensorgraph.models.graph_models import GraphConvModel
    model = GraphConvModel(n_tasks=1, batch_size=64, uncertainty=False, mode='regression')

    model = dc.models.GraphConvTensorGraph(
        1,
        batch_size=batch_size,
        learning_rate=1e-3,
        use_queue=False,
        mode = 'regression',
        model_dir= t_savename.get())

    np.random.seed(1)
    random.seed(1)
    
    model.fit(dataset, nb_epoch=max(1, int(t_epochs.get())))
    #model.fit(trainset, nb_epoch=max(1, int(t_epochs.get())))
    
    metric = dc.metrics.Metric(dc.metrics.r2_score)

    print('epoch: ', t_epochs.get() )
    print("Evaluating model")
    train_score = model.evaluate(trainset, [metric])
    test_score  = model.evaluate(testset, [metric])

    model.save()

    pred_train = model.predict(trainset)
    pred_test  = model.predict(testset)

    y_train = np.array(trainset.y, dtype = np.float32)
    y_test = np.array(testset.y, dtype = np.float32)
    
    import matplotlib.pyplot as plt

    plt.figure()

    plt.figure(figsize=(5,5))
    
    plt.scatter(y_train, pred_train, label = 'Train', c = 'blue')
    plt.title('Graph Convolution')
    plt.xlabel('Measured value')
    plt.ylabel('Predicted value')
    plt.scatter(y_test, pred_test, c = 'lightgreen', label = 'Test', alpha = 0.8)
    plt.legend(loc = 4)
    #plt.show()
    plt.savefig('score-tmp.png')
    
    from PIL import Image
    img = Image.open('score-tmp.png')

    img_resize = img.resize((400, 400), Image.LANCZOS)
    img_resize.save('score-tmp.png')


    global image_score
    image_score_open = Image.open('score-tmp.png')
    image_score = ImageTk.PhotoImage(image_score_open, master=frame1)

    canvas.create_image(200,200, image=image_score)

    #Calculate R2 score
    print("Train score")
    print(train_score)
    t_train_r2.set(train_score)

    print("Test scores")
    print(test_score)
    t_test_r2.set(test_score)

    #Calculate RMSE
    train_rmse = 1
    test_rmse  = 1
    '''
    print("Train RMSE")
    print(train_rmse)
    t_train_rmse.set(train_rmse)

    print("Test RMSE")
    print(test_rmse)
    t_test_rmse.set(test_rmse)
    '''

    df_save = pd.DataFrame(
    {'pred_train':pred_train,
    'meas_train':y_train
    })
    
    df_save.to_csv('pred_and_meas_train.csv')

    print('finish!')



def mordred_training():

    df_mordred = pd.read_csv(t_csv.get())

    def apply_molfromsmiles(smiles_name):
        try:
            mols = Chem.MolFromSmiles(smiles_name)
            
        except:
            mols = ""
            print(smiles_name)
            print('Error')

        return mols


    df_mordred['mols'] = df_mordred['smiles'].map(apply_molfromsmiles)
    df_mordred = df_mordred.dropna()


    calc =Calculator(descriptors, ignore_3D = True)            
    X = calc.pandas(df_mordred['mols']).astype('float').dropna(axis = 1)
    X = np.array(X, dtype = np.float32)

    #各記述子について平均0, 分散1に変換
    st = StandardScaler()
    X= st.fit_transform(X)
    
    #後で再利用するためにファイルに保存
    np.save("X_2d.npy", X)

    Y = df_mordred[t_task.get()]
    
    #Numpy形式の配列に変換
    Y = np.array(Y, dtype = np.float32)
    
    #後で再利用するためにファイルに保存
    np.save("Y_2d.npy", Y)

    X_train, X_test, y_train, y_test = model_selection.train_test_split(X,
    Y, test_size=0.25, random_state=42)

    model = Sequential()
    #入力層．Denseは全結合層の意味．次の層に渡される次元は50．入力データの次元（input_dim）は1114．
    model.add(Dense(units = 2400, input_dim = X.shape[1]))
    model.add(Activation("relu"))
    model.add(BatchNormalization())

    #出力層．次元1，つまり一つの値を出力する．
    model.add(Dense(units = 1))    
    model.summary()

    #SGDはStochastic Gradient Descent (確率的勾配降下法)．局所的最小値にとどまらないようにする方法らしい．nesterovはNesterovの加速勾配降下法．
    model.compile(loss = 'mean_squared_error',
    optimizer = 'Adam',
    metrics=['accuracy'])
    history = model.fit(X_train, y_train, epochs = max(5, int(t_epochs.get())), batch_size = 32,
    validation_data = (X_test, y_test))
    score = model.evaluate(X_test, y_test, verbose = 0)
    print('Test loss:', score[0])
    print('Test accuracy:', score[1])
    y_pred = model.predict(X_test)
    rms = (np.mean((y_test - y_pred) ** 2)) ** 0.5
    #s = np.std(y_test - y_pred)
    print("Neural Network RMS", rms)

    model.save('model_sol.h5')


    import matplotlib.pyplot as plt

    pred_train = model.predict(X_train, batch_size = 32)
    pred_test  = model.predict(X_test, batch_size = 32)


    df_save = pd.DataFrame(
        {'pred_train':pred_train,
        'meas_train':y_train
        })
    
    df_save.to_csv('pred_and_meas_train.csv')

    plt.figure(figsize=(5,5))
    plt.scatter(y_train, pred_train, label = 'Train', c = 'blue')
    plt.title('Mordred predict')
    plt.xlabel('Measured value')
    plt.ylabel('Predicted value')
    plt.scatter(y_test, pred_test, c = 'lightgreen', label = 'Test', alpha = 0.8)
    plt.legend(loc = 4)
    #plt.show()
    
    plt.savefig('score-tmp.png')

    global image_score
    image_score_open = Image.open('score-tmp.png')
    image_score = ImageTk.PhotoImage(image_score_open, master=frame1)

    canvas.create_image(200,200, image=image_score)


root = tkinter.Tk()

font1 = font.Font(family='游ゴシック', size=10, weight='bold')
root.option_add("*Font", font1)
root.option_add("*Button.font", font1)

style1  = ttk.Style()
style1.configure('my.TButton', font = ('游ゴシック',10)  )

root.title('test-horie')


frame1  = tkinter.ttk.Frame(root)

label_csv  = tkinter.ttk.Label(frame1, text = 'csv file path:')
t_csv = tkinter.StringVar()
entry_csv  = ttk.Entry(frame1, textvariable=t_csv, width = 60) 

button_getcsv = ttk.Button(frame1, text='CSVデータの選択', command = get_csv, style = 'my.TButton')

label_task = tkinter.ttk.Label(frame1, text = '特性値の列名:')
label_smiles = tkinter.ttk.Label(frame1, text = 'smilesの列名:')
label_id = tkinter.ttk.Label(frame1, text = 'IDの列名:')
label_savename = tkinter.ttk.Label(frame1, text = '保存フォルダ名:')

t_task = tkinter.StringVar()
t_smiles = tkinter.StringVar()
t_id = tkinter.StringVar()
t_savename = tkinter.StringVar()

entry_task = ttk.Entry(frame1, textvariable=t_task, width = 60) 
entry_smiles = ttk.Entry(frame1, textvariable=t_smiles, width = 60) 
entry_id = ttk.Entry(frame1, textvariable=t_id, width = 60) 
entry_savename = ttk.Entry(frame1, textvariable=t_savename, width = 60) 

t_epochs = tkinter.StringVar()
label_epochs = tkinter.ttk.Label(frame1, text = '学習回数:')
entry_epochs = ttk.Entry(frame1, textvariable=t_epochs, width = 60) 


button_graph_conv_training = ttk.Button(frame1, text = '畳み込みディープラーニング', command = graph_conv_training, style = 'my.TButton')
button_mordred_training = ttk.Button(frame1, text = 'mordred記述子　ディープラーニング', command = mordred_training, style = 'my.TButton')


label_train_r2 = tkinter.ttk.Label(frame1, text = '訓練用データのR2 score:')
t_train_r2 = tkinter.StringVar()
entry_train_r2 = ttk.Entry(frame1, textvariable=t_train_r2, width = 60) 

label_test_r2 = tkinter.ttk.Label(frame1, text = 'テスト用データのR2 score:')
t_test_r2 = tkinter.StringVar()
entry_test_r2 = ttk.Entry(frame1, textvariable=t_test_r2, width = 60) 

'''
label_train_rmse = tkinter.ttk.Label(frame1, text = '訓練用データのRMSE score:')
t_train_rmse = tkinter.StringVar()
entry_train_rmse = ttk.Entry(frame1, textvariable=t_train_rmse, width = 60) 

label_test_rmse = tkinter.ttk.Label(frame1, text = 'テスト用データのRMSE score:')
t_test_rmse = tkinter.StringVar()
entry_test_rmse = ttk.Entry(frame1, textvariable=t_test_rmse, width = 60) 
'''

frame1.grid(row=0,column=0,sticky=(N,E,S,W))
label_csv.grid(row=1,column=1,sticky=E)
entry_csv.grid(row=1,column=2,sticky=W)
button_getcsv.grid(row=2,column=2,sticky=W)

label_task.grid(row=4,column=1,sticky=E)
label_smiles.grid(row=5,column=1,sticky=E)
label_id.grid(row=6,column=1,sticky=E)
label_savename.grid(row=7,column=1,sticky=E)
label_epochs.grid(row=8,column=1,sticky=E)


entry_task.grid(row=4,column=2,sticky=W)
entry_smiles.grid(row=5,column=2,sticky=W)
entry_id.grid(row=6,column=2,sticky=W)
entry_savename.grid(row=7,column=2,sticky=W)
entry_epochs.grid(row=8,column=2,sticky=W)



button_graph_conv_training.grid(row=9,column=2,sticky=W)
button_mordred_training.grid(row=10,column=2,sticky=W)



label_train_r2.grid(row=11,column=1,sticky=E)
entry_train_r2.grid(row=11,column=2,sticky=W)

label_test_r2.grid(row=12,column=1,sticky=E)
entry_test_r2.grid(row=12,column=2,sticky=W)

#label_train_rmse.grid(row=13,column=1,sticky=E)
#entry_train_rmse.grid(row=13,column=2,sticky=W)

#label_test_rmse.grid(row=14,column=1,sticky=E)
#entry_test_rmse.grid(row=14,column=2,sticky=W)



canvas = tkinter.Canvas(frame1, width = 400, height = 400)
canvas.grid(row=13, column = 2)

global image_score

from PIL import ImageTk, Image
image_score_open = Image.open('deepchem_.png')

image_score = ImageTk.PhotoImage(image_score_open, master=frame1)
#image_ = image_.subsample(4)

canvas.create_image(200,200, image=image_score)


for child in frame1.winfo_children():
    child.grid_configure(padx=5, pady=5)

root.mainloop()