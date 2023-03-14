import numpy as np
from math import *
from interpolation import InterpVec

class Target(object):
       
    @classmethod
    def get_simple_target(cls, pos, vel):
        velocity_vectors = [[0, np.array(vel)]]
        vel_interp = InterpVec(velocity_vectors)
        target = cls(vel_interp=vel_interp)
        parameters_of_target = np.array([pos[0], pos[1], 0])
        target.set_init_cond(parameters_of_target=parameters_of_target)
        return target

    def __init__(self, *args, **kwargs):
        self.g   = kwargs.get('g', 9.80665)
        self.dt  = kwargs.get('dt', 0.001)
        self.vel_interp = kwargs['vel_interp']

    def set_init_cond(self, parameters_of_target=None):
        if parameters_of_target is None:
            parameters_of_target = self.get_standart_parameters_of_target()
        self.state = np.array(parameters_of_target)
        self.state_0 = np.array(parameters_of_target)

    def reset(self):
        self.set_state(self.state_0)

    def set_state(self, state):
        self.state = np.array(state)

    def get_state(self):
        return self.state
    
    def get_state_0(self):
        return self.state_0

    def step(self, tau):
        x, y, t = self.state
        t_end = t + tau
        flag = True
        while flag:
            if t_end - t > self.dt:
                dt = self.dt 
            else:
                dt = t_end - t
                flag = False
            t += dt
            vx, vy = self.vel_interp(t)
            x += vx * dt 
            y += vy * dt
        self.set_state([x, y, t])

    @property
    def pos(self):
        return self.state[:2]
    
    @property
    def vel(self):
        return self.vel_interp(self.t)

    @property
    def t(self):
        return self.state[-1]
    
    @property
    def Q(self):
        vx, vy = self.vel_interp(self.t)
        return np.sqrt(vx ** 2 + vy ** 2)

    @property
    def v(self):
        vx, vy = self.vel_interp(self.t)
        return np.sqrt(vx ** 2 + vy ** 2)

    @property
    def x(self):
        return self.pos[0]

    @property
    def y(self):
        return self.pos[1]

    def get_summary(self):
        return { 
            't': self.t,
            'v': self.v,
            'x': self.x,
            'y': self.y,
            'Q': np.degrees(self.Q)
        }
    
import numpy as np
from easyvec import Vec2


class Target2D(object):

    @classmethod
    def get_simple_target(cls, opts, time_min, postProcessing, g=9.80665, dt=0.01):
        pos, pos_aim, vel, vel_aim = opts['pos'], opts['pos_aim'], opts['vel'], opts['vel_aim'],
        pos = cls._validate_dimension_argument(pos)
        vel = cls._validate_dimension_argument(vel)
        aim = cls._validate_dimension_argument(pos_aim)
        vel_aim = cls._validate_dimension_argument(vel_aim)
        target = cls(pos=pos, vel=vel, aim=aim, vel_aim=vel_aim, g=g, dt=dt,
                     time_min=time_min, postProcessing=postProcessing)
        parameters_simple_target = np.array([pos[0], pos[1], vel[0], vel[1], 0])
        target.set_initial_condition(init_parameters=parameters_simple_target)
        return target

    def __init__(self, **kwargs):
        self.g = kwargs['g']
        self.t = 0
        self.dt = kwargs['dt']
        
        self.A = kwargs['pos']
        self.D = kwargs['aim']
        self.velA = kwargs['vel'] 
        self.velD = kwargs['vel_aim']
        
        self.fly_time = 0
        self._set_fly_time(kwargs['time_min'])
        
        self.n = 5 
        self.state = np.zeros(self.n)
        self.state_init = np.zeros(self.n)
        
        self.postProcessing = kwargs['postProcessing']
        self.history = {
                'x': [],
                'y': [],
                'v_abs': [],
                'vx': [],
                'vy': [],
                't': []
            }
    
    @property
    def pos(self):
        return self._fpos(self.t)
    
    @property
    def x(self):
        return self._fpos(self.t)[0]

    @property
    def y(self):
        return self._fpos(self.t)[1]
    
    @property
    def vel(self):
        return self._fvel(self.t)
    
    @property
    def vel_abs(self):
        return np.sqrt(self.vel.dot(self.vel))
    
    @property
    def vel_x(self):
        return self.vel[0]
    
    @property
    def vel_y(self):
        return self.vel[1]
    
    def set_initial_condition(self, init_parameters=None):
        if init_parameters is None:
            init_parameters = self.get_random_parameters()
        self.state = np.array(init_parameters)
        self.state_init = np.array(init_parameters)

    def get_random_parameters(self):
        # TODO
        return 0
    
    def step(self, tau=0.1, n=10):
        
        t = self.t  
        dt = tau / n
        
        for i in range(n):
            t += dt
            self.t = t
            state = [self.x, self.y, self.vel_x, self.vel_y, t]
            self.set_state(state)
            
            if self.postProcessing:
                self.history['x'].append(self.x)
                self.history['y'].append(self.y)
                self.history['v_abs'].append(self.vel_abs)
                self.history['vx'].append(self.vel_x)
                self.history['vy'].append(self.vel_y)
                self.history['t'].append(self.t)
                
    def reset(self):
        self.set_state(self.state_init)

    def set_state(self, state):
        self.state = np.array(state)

    def get_state(self):
        return self.state
    
    def get_state_init(self):
        return self.state_init
    
    def get_B(self, fly_time):
        return self.A + (self.velA * fly_time / 3)
    
    def get_C(self, fly_time): 
        return self.D - (self.velD * fly_time / 3)
    
    def get_traject(self, fly_time, vel_aim=None, n_points=100):
        vel_aim = self.velD if vel_aim is None else vel_aim
        taus = np.linspace(0, 1, n_points)
        A, velA = self.A, self.velA
        D, velD = self.D, vel_aim
        B = self.get_B(fly_time)
        C = self.get_C(fly_time)
        return np.array([
            (1-tau)**3 * A + 3*tau*(1-tau)**2 * B + 3*tau*tau*(1-tau)*C + tau**3 * D
            for tau in taus])
    
    def get_traject_vels(self, fly_time, vel_aim=None, n_points=100):
        vel_aim = self.velD if vel_aim is None else vel_aim
        taus = np.linspace(0, 1, n_points)
        A, velA = self.A, self.velA
        D, velD = self.D, vel_aim
        B = self.get_B(fly_time)
        C = self.get_C(fly_time)
        return np.array([
            (3*(1-tau)**2*(B-A) + 6*tau*(1-tau)*(C-B) + 3*tau**2*(D-C)) / fly_time
            for tau in taus
        ])
    
    def get_traject_acc(self, fly_time, vel_aim=None, n_points=100):
        vel_aim = self.velD if vel_aim is None else vel_aim
        taus = np.linspace(0, 1, n_points)
        A, velA = self.A, self.velA
        D, velD = self.D, vel_aim
        B = self.get_B(fly_time)
        C = self.get_C(fly_time)
        return np.array([
            6 * ((tau-1)*(B-A) + (1-2*tau)*(C-B) + tau*(D-C)) / (fly_time**2)
            for tau in taus
        ])
    
    def get_fly_time_minimum(self, ny):
        # TODO: помаксимальной перегрузке вычислить время
        # n = overload if overload != None else self.overload
        #         for i in :
        #             a_max = self.get_acc_max(fly_time)
        #             if (a_max / self.g) > n_max:
        #                 n_max = a_max / self.g
        #         return
        return 0
    
    def get_fly_time_random(self):
        fly_time_A = ((self.A - self.D).len() / self.velA.len())
        fly_time_D = (self.A - self.D).len() / self.velD.len()
        fly_time_min = min(fly_time_A, fly_time_D)
        fly_time_max = max(fly_time_A, fly_time_D)
        return np.random.uniform(fly_time_min, fly_time_max)

    def get_vel_min(self, fly_time, n=51):
        BA = self.B.sub_vec(self.A)
        CB = self.C.sub_vec(self.B)
        DC = self.D.sub_vec(self.C)
        min_v = 1e6
        dt = 1.0 / (n - 1)
        for i in range(n):
            t = i * dt
            vel_len = (BA.mul_num(3*(1-t)*(1-t)/fly_time).add_vec(CB.mul_num(6*t*(1-t)/fly_time)).add_vec(DC.mul_num(3*t*t/fly_time))).len()
            if vel_len < min_v:
                min_v = vel_len
        return min_v    
    
    def get_vel_max(self, delta_t, n=51):
        BA = self.B.sub_vec(self.A)
        CB = self.C.sub_vec(self.B)
        DC = self.D.sub_vec(self.C)
        max_v = 0
        dt = 1.0 / (n - 1)
        for i in range(n):
            t = i * dt
            vel_len = (BA.mul_num(3*(1-t)*(1-t)/delta_t).add_vec(CB.mul_num(6*t*(1-t)/delta_t)).add_vec(DC.mul_num(3*t*t/delta_t))).len() 
            if vel_len > max_v:
                max_v = vel_len
        return max_v   
    
    def get_acc_max(self, fly_time=None):
        fly_time = fly_time if fly_time is not None else self.fly_time
        B, C = self.get_B(fly_time), self.get_C(fly_time)
        a1 = (C - B) * 6 - (B - self.A) * 6 - self.G
        a2 = (self.D - C) * 6 - (C - B) * 6 - self.G
        return np.fmax(a1.len(), a2.len()) / (fly_time**2)
    
    def _fpos(self, t, fly_time=None):
        fly_time = fly_time if fly_time is not None else self.fly_time
        A, velA = self.A, self.velA
        D, velD = self.D, self.velD
        B = self.get_B(fly_time)
        C = self.get_C(fly_time)
        tau = t / fly_time
        return (1-tau)**3 * A + 3*tau*(1-tau)**2 * B + 3*tau*tau*(1-tau)*C + tau**3 * D
    
    def _fvel(self, t, fly_time=None):
        fly_time = fly_time if fly_time is not None else self.fly_time
        A, velA = self.A, self.velA
        D, velD = self.D, self.velD
        B = self.get_B(fly_time)
        C = self.get_C(fly_time)
        tau = t / fly_time
        return (3*(1-tau)**2*(B-A) + 6*tau*(1-tau)*(C-B) + 3*tau**2*(D-C)) / fly_time
    
    def _facc(self, t, fly_time=None):
        fly_time = fly_time if fly_time is not None else self.fly_time
        A, velA = self.A, self.velA
        D, velD = self.D, self.velD
        B = self.get_B(fly_time)
        C = self.get_C(fly_time)
        tau = t / self.fly_time 
        return 6 * ((tau-1)*(B-A) + (1-2*tau)*(C-B) + tau*(D-C)) / (fly_time**2)
 
    def to_dict(self):
        return { 
            't': self.t,
            'x': self.x,
            'y': self.y,
            'v': self.vel,
        }
    
    def _set_fly_time(self, time_min=True):
        self.fly_time = self.get_fly_time_minimum() if time_min else self.get_fly_time_random()
    
    @staticmethod
    def _validate_dimension_argument(array, n=2):
        if len(array) == 0:
            return Vec2(0, 0)
        elif len(array) == n:
            try:
                l = [float(elem) for elem in array]
            except ValueError:
                raise ValueError("Один или несколько элементов в <{!r}> не действительное(-ые) число(-а)".format(seq))
            else:
                arg1, arg2 = l
                return Vec2(arg1, arg2)  
        else:
            raise ValueError("Неожиданное число элементов. Получено: {}, Ожидалось: {}.".format(len(array), n))  