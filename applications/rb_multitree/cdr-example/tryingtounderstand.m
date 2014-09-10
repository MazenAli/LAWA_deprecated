format long



load('Runs/Stage7_Run3/training_data_stage7/RB_A_1.txt');
load('Runs/Stage7_Run3/training_data_stage7/RB_A_2.txt');
load('Runs/Stage7_Run3/training_data_stage7/RB_A_3.txt');
load('Runs/Stage7_Run3/training_data_stage7/RB_A_4.txt');
load('Runs/Stage7_Run3/training_data_stage7/RB_F_1.txt');

BF1 = importdata('Runs/Stage7_Run3/training_data_stage7/bf/bf_1.txt', ' ', 1);
BF2 = importdata('Runs/Stage7_Run3/training_data_stage7/bf/bf_2.txt', ' ', 1);
BF3 = importdata('Runs/Stage7_Run3/training_data_stage7/bf/bf_3.txt', ' ', 1);
BF4 = importdata('Runs/Stage7_Run3/training_data_stage7/bf/bf_4.txt', ' ', 1);
BF5 = importdata('Runs/Stage7_Run3/training_data_stage7/bf/bf_5.txt', ' ', 1);
BF6 = importdata('Runs/Stage7_Run3/training_data_stage7/bf/bf_6.txt', ' ', 1);

BF1 = BF1.data;
BF2 = BF2.data;
BF3 = BF3.data;
BF4 = BF4.data;
BF5 = BF5.data;
BF6 = BF6.data;


A = RB_A_1 + RB_A_2 +    3.6936058457034515e-01 * RB_A_3 + 3.1396037564550863 * RB_A_4
F = RB_F_1'
 
A4 = A(1:4, 1:4);
F4 = F(1:4);
u4 = A4 \ F4
A4*u4 - F4

A5 = A(1:5, 1:5);
F5 = F(1:5);
u5 = A5 \ F5
A5*u5 - F5

upr_4 = A(2:5, 2:5)\F(2:5)
A(2:5, 2:5) * upr_4 - F(2:5)

u5_pr = [0;  upr_4] 
A5 * u5_pr - F5

