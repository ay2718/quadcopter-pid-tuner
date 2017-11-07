import math
import numpy as np
import matplotlib.pyplot as plt

stab_p = 0.0;
stab_d = 0.0;
KP = 0.05;
KI = 0.2;
KD = 0.002;
bias = 0.0;
noise = False;

yn = input("Would you like to use default values (Y/n)? ")
if yn == 'n' or yn == 'N':
    KP = float(input("KP = "));
    KI = float(input("KI = "));
    KD = float(input("KD = "));
else:
    print("KP =", KP);
    print("KI =", KI);
    print("KD =", KD);

yn = input("Would you like to stabilize position (y/N)? ")
if yn == 'y' or yn == 'Y':
    stab_p = float(input("Position KP = "));
#     stab_d = float(input("Position KD = "));

yn = input("Would you like to add a systematic error (y/N)? ")
if yn == 'y' or yn == 'Y':
    bias = 0.02;

yn = input("Would you like to add random noise (y/N)? ")
if yn == 'y' or yn == 'N':
    noise = True;

base_thrust = 0.5;

time = 5.0
voltage = 10.0
motor_kv = 1090.0

motor_kv *= (math.pi/30.0)
motor_r = 0.235
quad_moment = 0.004
prop_moment = 5e-5
arm_length = 0.15

quad_speed = 0.0
quad_position = 0.0

# step size of physics
h = 0.001

steps = int(time//h);
quad_speed_log = np.empty(steps);
quad_position_log = np.empty(steps);
prop_speed_log = np.empty((steps, 2));
prop_voltage_log = np.empty((steps, 2));
setpoint_log = np.empty(steps);

# step size of control circuitry
arduino_h = 0.01
arduino_steps = int(arduino_h // h);

# delay between changing motor speed and receiving data
delay = 0.03
delay_steps = int(delay // h);

prop_speed = np.zeros(2)
prop_voltage = np.zeros(2)

# prop_thr_const satisfies thrust = w^2 * prop_thr_const
prop_thr_const = 3e-7
# prop_res_const satisfies torque = w^2 * prop_res_const
prop_res_const = 1e-7

thrust = np.zeros(2);

def thrust_to_voltage(): return voltage * np.sqrt(thrust)

# gives propeller acceleration as a function of a bunch of parameters
def prop_acc():
#    we know that I = (V - w/Kv)/R and that total power = IV = V(V - w/Kv)/R
#    usable power is w/Kv*(V - w/Kv)/R = 
#    usable power = torque*w, so torque = (V - w/Kv)/(R*Kv)
#    counter torque from the prop in w^2 * rconst
#    Therefore, total torque is (V - w/Kv)/(R*Kv) - rconst*w^2
#    a = torque/moment
#    Also note that the motor freewheels when w/Kv > V, so we have the following formula
    global prop_voltage;
    pac = np.zeros(2);
    for i in range(0, 2):
        if prop_speed[i] > prop_voltage[i]*motor_kv:
            pac[i] = -prop_res_const*prop_speed[i]**2/prop_moment
        else:
            pac[i] = ((prop_voltage[i] - prop_speed[i]/motor_kv)/(motor_kv*motor_r) - prop_res_const*prop_speed[i]**2)/prop_moment
    return pac;

def quad_acc():
    return 2*(prop_speed[1]**2 - prop_speed[0]**2)*prop_thr_const*arm_length/prop_moment

def physics_step(i):
    global quad_speed
    global quad_position
    global prop_speed
    quad_speed_log[i] = quad_speed;
    quad_position_log[i] = quad_position;
    prop_speed_log[i, :] = prop_speed;
    prop_voltage_log[i, :] = prop_voltage;
    setpoint_log[i] = setpoint(i);

    kprop = np.empty((4, 2));
    kquadspeed = np.empty(4);
    kquadpos = np.empty(4);

    kprop[0, :] = prop_acc();
    kquadspeed[0] = quad_acc();
    kquadpos[0] = quad_speed;

    prop_speed = prop_speed_log[i, :] + kprop[0, :]*h/2;
    quad_speed = quad_speed_log[i] + kquadspeed[0]*h/2;
    quad_position = quad_position_log[i] + kquadpos[0]*h/2;

    kprop[1, :] = prop_acc();
    kquadspeed[1] = quad_acc();
    kquadpos[1] = quad_speed;

    prop_speed = prop_speed_log[i, :] + kprop[1, :]*h/2;
    quad_speed = quad_speed_log[i] + kquadspeed[1]*h/2;
    quad_position = quad_position_log[i] + kquadpos[1]*h/2;

    
    kprop[2, :] = prop_acc();
    kquadspeed[2] = quad_acc();
    kquadpos[2] = quad_speed;

    prop_speed = prop_speed_log[i, :] + kprop[2, :]*h;
    quad_speed = quad_speed_log[i] + kquadspeed[2]*h;
    quad_position = quad_position_log[i] + kquadpos[2]*h;


    kprop[3, :] = prop_acc();
    kquadspeed[3] = quad_acc();
    kquadpos[3] = quad_speed;

    prop_speed = prop_speed_log[i, :] + (kprop[0, :] + kprop[3, :] + 2*(kprop[1, :] + kprop[2, :]))*h/6;
    quad_speed = quad_speed_log[i] + (kquadspeed[0] + kquadspeed[3] + 2*(kquadspeed[1] + kquadspeed[2]))*h/6;
    quad_position = quad_position_log[i] + (kquadpos[0] + kquadpos[3] + 2*(kquadpos[1] + kquadpos[2]))*h/6;

shelf = 2.0
def setpoint(curr_step):
    if curr_step*h > shelf:
        return 1.0;
    else:
        return 0.0;

integral_error = 0.0;

def error(i):
    global integral_error;
    qspeed = quad_speed;
    if noise: qspeed += np.random.normal(0.0, 0.001);
    err = 0.0;
    if stab_p != 0:
        err = qspeed - stab_p*(setpoint(i) - quad_position) + stab_d*quad_speed;
    else:
        err = qspeed - setpoint(i);
    k = i - arduino_steps
    if k < 0: k = 0;
    integral_error += arduino_h*err
    derivative_error = (qspeed - quad_speed_log[k])/arduino_h;
    correction = KP*err + KI*integral_error + KD*derivative_error;
    return correction;
    
def trim(value):
    if value > 1.0: value = 1.0;
    if value < 0.0: value = 0.0;
    return value;

for i in range(0, steps):
    physics_step(i);
    if i % arduino_steps == 0:
        j = i - delay_steps;
        if j < 0: j = 0;
        correction = error(i) + bias;
        thrust[0] = base_thrust + correction;
        thrust[1] = base_thrust - correction;
        thrust[0] = trim(thrust[0]);
        thrust[1] = trim(thrust[1]);
        prop_voltage = thrust_to_voltage();

_, plots = plt.subplots(2, sharex = True)
plots[0].plot(h*np.arange(0, steps), quad_speed_log);
if stab_p == 0:
    plots[0].plot(h*np.arange(0, steps), setpoint_log);

plots[1].plot(h*np.arange(0, steps), quad_position_log);
if stab_p != 0:
    plots[1].plot(h*np.arange(0, steps), setpoint_log);
plt.show();
