import math
import numpy as np
import matplotlib.pyplot as plt

KP = 0.1;
KI = 0.5;
KD = 0.0001;

ynstr = input("Would you like to tune your motor? (n --> no tuning, p --> P only, i --> PI tuning) ");
if ynstr == 'n' or ynstr == 'N':
    KP = KI = KD = 0.0;
elif ynstr == 'p' or ynstr == 'P':
    KI = KD = 0.0;
else:
    KD = 0.0;
# yn = input("Would you like to use default values (Y/n)? ")
# if yn == 'n' or yn == 'N':
#     KP = float(input("KP = "));
#     KI = float(input("KI = "));
#     KD = float(input("KD = "));
# else:
#     print("KP =", KP);
#     print("KI =", KI);
#     print("KD =", KD);

time = 5.0
voltage = 10.0
motor_kv = 1090.0

motor_kv *= (math.pi/30.0)
motor_r = 0.235
motor_moment = 1e-4

# step size of physics
h = 0.001

steps = int(time//h);
motor_speed_log = np.empty(steps);
motor_voltage_log = np.empty(steps);
setpoint_log = np.empty(steps);
resist_log = np.empty(steps);

# step size of control circuitry
arduino_h = 0.01
arduino_steps = int(arduino_h // h);

# delay between changing motor speed and receiving data
delay = 0.03
delay_steps = int(delay // h);

motor_speed = 0.0
motor_voltage = 0.0

# gives motor acceleration as a function of a bunch of parameters
def motor_acc(resistance):
    global motor_voltage;
    global motor_speed;
    acc = 0.0
    if motor_speed > motor_voltage*motor_kv:
        acc = -resistance/motor_moment
    else:
        acc = ((motor_voltage - motor_speed/motor_kv)/(motor_kv*motor_r) - resistance)/motor_moment
    if acc < 0 and motor_speed <= 0:
        return 0
    else:
        return acc

def quad_acc():
    return 2*(prop_speed[1]**2 - prop_speed[0]**2)*prop_thr_const*arm_length/prop_moment

def physics_step(i):
    global motor_speed
    motor_speed_log[i] = motor_speed;
    motor_voltage_log[i] = motor_voltage;
    setpoint_log[i] = setpoint(i);
    resist_log[i] = resist(i);

    kmotor = np.empty(4);

    kmotor[0] = motor_acc(resist_log[i]);
    motor_speed = motor_speed_log[i] + kmotor[0]*h/2;
    
    kmotor[1] = motor_acc(resist_log[i]);
    motor_speed = motor_speed_log[i] + kmotor[1]*h/2;
    
    kmotor[2] = motor_acc(resist_log[i]);
    motor_speed = motor_speed_log[i] + kmotor[2]*h;

    kmotor[3] = motor_acc(resist_log[i]);
    motor_speed = motor_speed_log[i] + (kmotor[0] + kmotor[3] + 2*(kmotor[1] + kmotor[2]))*h/6;

shelf = 0.5
def setpoint(curr_step):
    setp = 0.0;
    if curr_step*h > shelf:
        setp = 1.0;
    else:
        setp = 0.0;
    return 600*setp;

shelf2 = 3.0
def resist(curr_step):
    res = 0.0
    if curr_step*h > shelf2:
        res = 3.0;
    else:
        res = 1.5;
    return 0.03*res;

integral_error = 0.0;

def error(i):
    global integral_error;
    err = setpoint(i) - motor_speed;
    k = i - arduino_steps
    if k < 0: k = 0;
    if err < 100: integral_error += arduino_h*err
    derivative_error = (motor_speed_log[k] - motor_speed)/arduino_h;
    correction = KP*err + KI*integral_error + KD*derivative_error;
    return correction;
    
def trim(value):
    if value > 10.0: value = 10.0;
    if value < 0.0: value = 0.0;
    return value;

for i in range(0, steps):
    physics_step(i);
    if i % arduino_steps == 0:
        j = i - delay_steps;
        if j < 0: j = 0;
        motor_voltage = trim(setpoint(i)/motor_kv + error(i));

_, plots = plt.subplots(2, sharex = True);
plots[0].plot(h*np.arange(0, steps), motor_speed_log);
plots[0].plot(h*np.arange(0, steps), setpoint_log);
plots[1].plot(h*np.arange(0, steps), motor_voltage_log);
plt.show();
