import math
import numpy as np
import matplotlib.pyplot as plt

KP = 40.0;
KI = 200.0;
KD = 4.0;

try: KP = float(input("KP = (default is "+str(KP)+") "));
except: pass
try: KI = float(input("KP = (default is "+str(KI)+") "));
except: pass
try: KD = float(input("KP = (default is "+str(KD)+") "));
except: pass

time = 5.0
voltage = 10.0
motor_kv = 30.0

motor_kv *= (math.pi/30.0)
motor_r = 0.235
arm_length = 0.25
arm_sym_mass = 2.0
arm_assym_mass = 0.0
try: arm_assym_mass = float(input("Enter a value for assymetrical arm mass in kg (default is 0) "));
except: pass
arm_mass = abs(arm_assym_mass) + arm_sym_mass
motor_moment = arm_length**2*arm_mass

# step size of physics
h = 0.001

steps = int(time//h);
motor_speed_log = np.empty(steps);
motor_voltage_log = np.empty(steps);
motor_position_log = np.empty(steps);
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
motor_position = 0.0

# gives motor acceleration as a function of a bunch of parameters
def motor_acc(resistance):
    global motor_voltage;
    global motor_speed;
    return ((motor_voltage - motor_speed/motor_kv)/(motor_kv*motor_r) - resistance)/motor_moment

def physics_step(i):
    global motor_speed
    global motor_position
    motor_speed_log[i] = motor_speed;
    motor_position_log[i] = motor_position;
    motor_voltage_log[i] = motor_voltage;
    setpoint_log[i] = setpoint(i);
    resist_log[i] = 10*arm_length*arm_assym_mass*math.sin(motor_position);

    kmotor = np.empty(4);
    kmotorpos = np.empty(4);

    kmotor[0] = motor_acc(resist_log[i]);
    kmotorpos[0] = motor_speed;
    motor_speed = motor_speed_log[i] + kmotor[0]*h/2;
    motor_position = motor_position_log[i] + kmotorpos[0]*h/2;
    
    kmotor[1] = motor_acc(resist_log[i]);
    kmotorpos[1] = motor_speed;
    motor_speed = motor_speed_log[i] + kmotor[1]*h/2;
    motor_position = motor_position_log[i] + kmotorpos[1]*h/2;
    
    kmotor[2] = motor_acc(resist_log[i]);
    kmotorpos[2] = motor_speed;
    motor_speed = motor_speed_log[i] + kmotor[2]*h;
    motor_position = motor_position_log[i] + kmotorpos[2]*h/2;

    kmotor[3] = motor_acc(resist_log[i]);
    kmotorpos[3] = motor_speed;
    motor_speed = motor_speed_log[i] + (kmotor[0] + kmotor[3] + 2*(kmotor[1] + kmotor[2]))*h/6;
    motor_position = motor_position_log[i] + (kmotorpos[0] + kmotorpos[3] + 2*(kmotorpos[1] + kmotorpos[2]))*h/6;

shelf = 1.0 
def setpoint(curr_step):
    setp = 0.0;
    if curr_step*h > shelf:
        setp = math.pi/2;
    else:
        setp = 0.0;
    return setp;

integral_error = 0.0;

def error(i):
    global integral_error;
    err = setpoint(i) - motor_position;
    k = i - arduino_steps
    if k < 0: k = 0;
    derivative_error = (motor_position_log[k] - motor_position)/arduino_h;
    if abs(KD*derivative_error) < 10.0 and abs(KP*err) < 10.0: integral_error += arduino_h*err
    correction = KP*err + KI*integral_error + KD*derivative_error;
    return correction;
    
def trim(value):
    if value > 10.0: value = 10.0;
    if value < -10.0: value = -10.0;
    return value;

for i in range(0, steps):
    physics_step(i);
    if i % arduino_steps == 0:
        j = i - delay_steps;
        if j < 0: j = 0;
        motor_voltage = error(i);
        motor_voltage = trim(motor_voltage);

_, plots = plt.subplots(2, sharex = True);
plots[0].plot(h*np.arange(0, steps), motor_position_log);
plots[0].plot(h*np.arange(0, steps), setpoint_log);
plots[1].plot(h*np.arange(0, steps), motor_voltage_log);
plots[0].set_title("Motor Position vs Time")
plots[1].set_title("Motor Voltage vs Time")
plt.show();
