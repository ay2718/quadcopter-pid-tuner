#
#
#
#
# Author: Aaron Yeiser
# PID.py simulates a quadcopter that is being PID tuned
# The user can also add noise or systematic error
#
#
# 

import math
import numpy as np
import matplotlib.pyplot as plt

# PID coefficients for position
stab_p = 0.0;
stab_d = 0.0;

# PID coefficients for angular velocity
KP = 0.05;
KI = 0.2;
KD = 0.002;

# bias is systematic error, noise is noise level in gyroscope
bias = 0.0;
noise = False;

# Allow user to change PID coefficients from defaults
try: KP = float(input("KP = (default is "+str(KP)+") "));
except: pass
try: KI = float(input("KP = (default is "+str(KI)+") "));
except: pass
try: KD = float(input("KP = (default is "+str(KD)+") "));
except: pass

# User can add systematic error (default is no error)
yn = input("Would you like to add a systematic error (y/N)? ")
if yn == 'y' or yn == 'Y':
    bias = 0.02;

# User can add random noise (default is no noise)
yn = input("Would you like to add random noise (y/N)? ")
if yn == 'y' or yn == 'N':
    noise = True;

# User can force quadcopter to add an addition position stabilization loop (default is no)
try: stab_p = float(input("Position stabilization P coefficient (default is 0) "));
except: pass

# Motors operate at 50% thrust by default
base_thrust = 0.5;

# Simulation runs for time, motors operate at voltage, and motor rotates at kv rpm per volt under no load
time = 5.0
voltage = 10.0
motor_kv = 1090.0

# convert to rad/s/volt
motor_kv *= (math.pi/30.0)

# motor resistance, moments of inertia
motor_r = 0.235
quad_moment = 0.004
prop_moment = 5e-5
arm_length = 0.15

quad_speed = 0.0
quad_position = 0.0

# step size of physics
h = 0.001

# Logging
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

# returns quadcopter angular acceleration
def quad_acc():
    return 2*(prop_speed[1]**2 - prop_speed[0]**2)*prop_thr_const*arm_length/prop_moment

# coupled Runge-Kutta solve for quadcopter position and speed and propeller speed
def physics_step(i):
    global quad_speed
    global quad_position
    global prop_speed

    # Log current prop speed values
    # Log values get used as temp variables here
    quad_speed_log[i] = quad_speed;
    quad_position_log[i] = quad_position;
    prop_speed_log[i, :] = prop_speed;
    prop_voltage_log[i, :] = prop_voltage;
    setpoint_log[i] = setpoint(i);

    # Runge Kutta magic happens
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

# setpoint function to force quadcopter to rotate at a specific time
shelf = 2.0
def setpoint(curr_step):
    if curr_step*h > shelf:
        return 1.0;
    else:
        return 0.0;

integral_error = 0.0;

# gives quadcopter KP*proportional error + KI*integral error + KD*derivative error
def error(i):
    global integral_error;
    qspeed = quad_speed;

    # Allows for noise to be simulated
    if noise: qspeed += np.random.normal(0.0, 0.001);
    err = 0.0;

    # Different methods of stabilizing the quadcopter
    if stab_p != 0:
        err = qspeed - stab_p*(setpoint(i) - quad_position) + stab_d*quad_speed;
    else:
        err = qspeed - setpoint(i);

    # Introduce delay
    k = i - arduino_steps
    if k < 0: k = 0;

    # compute error
    integral_error += arduino_h*err
    derivative_error = (qspeed - quad_speed_log[k])/arduino_h;
    correction = KP*err + KI*integral_error + KD*derivative_error;
    return correction;
    
# trims a value to between 0 and 1
def trim(value):
    if value > 1.0: value = 1.0;
    if value < 0.0: value = 0.0;
    return value;

# simulates quadcopter
for i in range(0, steps):
    physics_step(i);
    if i % arduino_steps == 0:
        j = i - delay_steps;
        if j < 0: j = 0;

        # add systematic error, if present
        correction = error(i) + bias;
        thrust[0] = base_thrust + correction;
        thrust[1] = base_thrust - correction;
        thrust[0] = trim(thrust[0]);
        thrust[1] = trim(thrust[1]);
        prop_voltage = thrust_to_voltage();

_, plots = plt.subplots(2, sharex = True)
plots[0].set_title("Quadcopter Angular Velocity vs Time");
plots[1].set_title("Quadcopter Angular Position vs Time");
plots[0].plot(h*np.arange(0, steps), quad_speed_log);
if stab_p == 0:
    plots[0].plot(h*np.arange(0, steps), setpoint_log);

plots[1].plot(h*np.arange(0, steps), quad_position_log);
if stab_p != 0:
    plots[1].plot(h*np.arange(0, steps), setpoint_log);
plt.show();
