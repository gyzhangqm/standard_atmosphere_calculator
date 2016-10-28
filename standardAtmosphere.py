#Atmospheric pressure and air density calculator

import math
import simplegui

#Boundary conditions at sea level

# Pressure and temperature at the sea level /assumed/
P0 = 101325.0
T0 = 288.15
rho0 = 1.225

g = 9.80665
# Ideal gas constant for air = ideal gas constant / molar mass of air
R = 287.0
X = 0.0
zone = 0
zonelist = [0]

# Lapse rate = temperature change in K per altitude change in km
ares = 0.0
a = 0.0
a1 = -6.5
a2 = 0.0
a3 = 1.0
a4 = 2.8
a5 = 0.0
a6 = -2.8
a7 = -2.0

# Zone altitude
hprev = 0.0
h0 = 0.0
h1 = 11.0
h2 = 20.0
h3 = 32.0
h4 = 47.0
h5 = 51.0
h6 = 71.0
h7 = 84.0

Tprev = T0
Pprev = P0
Tres = T0
Pres = P0
rhores = rho0
Tdifference = 0
Pper = 100.0
rhoper = 100.0


def draw(canvas):
    canvas.draw_polygon([(2,135), (2, 265), (400,265), (400,135)], 1, 'Navy')
    canvas.draw_text('This program calculates temperature, pressure and air density', (4,20), 15, 'Black')
    canvas.draw_text('at altitudes up to 84 km above the sea level.', (4, 40), 15, 'Black')
    if X <= 84 and X >=0 :
        canvas.draw_text('Your altitude is: ' + str(X) + ' km', (4,100), 15, 'Navy')
        canvas.draw_text('Temperature at altitude ' + str(X) + ' km: ' + str(round(Tres, 3)) + ' K.', (4,150), 15, 'Black')
        canvas.draw_text('Pressure at altitude ' + str(X) + ' km: ' + str(round(Pres,3)) + ' Pa.', (4,200), 15, 'Black')
        canvas.draw_text('Air density at altitude ' + str(X) + ' km: ' + str(round(rhores,3)) + ' kg/m3.', (4,245), 15, 'Black')
        
        canvas.draw_text('Temperature difference with respect to sea level: ' + str(round(Tdifference,2)) + ' deg.', (4,165) , 15, 'Black' )
        canvas.draw_text('Percentage of sea level pressure: ' + str(round(Pper,3)) + ' %.', (4,215),15 ,'Black')    
        canvas.draw_text('Percentage of sea level density: ' + str(round(rhoper,3)) + ' %.', (4,260), 15, 'Black')    

    if X == 'NaN':
        canvas.draw_text('Value entered is not a number.', (4,100), 15, 'Navy')
    if X == -1:
        canvas.draw_text('Altitude cannot be greater than 84 km or negative', (4, 120), 15, 'Black')
    
  
def inp(number):
    global X, zone, a, isothermal, zonelist, ares, hprev, Pres, rhores, Tdifference, Pper, rhoper, Tres
    
    try:
        float(number)
        check = True
        X = float(number)
    except ValueError:
        check = False
    
    if check == True and X <= 84 and X >=0:

        if X >= h0 and X <= h1:
            zone = 1
            a = a1
            zonelist = [1]
            
            Tprev = T0
            ares = a1
            hprev = h0
            
        elif X > h1 and X <= h2:
            zone = 2
            a = a2
            zonelist = [1,2]
            
            ares = a2
            hprev = h1
            
        elif X > h2 and X <= h3:
            zone = 3
            a = a3
            zonelist = [1,2,3]
            
            ares = a3
            hprev = h2
            
        elif X > h3 and X <= h4:
            zone = 4
            a = a4
            zonelist = [1,2,3,4]
            
            ares = a4
            hprev = h3
            
        elif X > h4 and X <= h5:
            zone = 5
            a = a5
            zonelist = [1,2,3,4,5]
            
            ares = a5
            hprev = h4
            
        elif X > h5 and X <= h6:
            zone = 6
            a = a6
            zonelist = [1,2,3,4,5,6]
            
            ares = a6
            hprev = h5
            
        elif X > h6 and X <= h7:
            zone = 7
            a = a7
            zonelist = [1,2,3,4,5,6,7]
            
            ares = a7
            hprev = h6
            
        else:
            zone = 0
            a = 0
            zonelist = [0]
    
        for i in zonelist:
            if i == 1:
                T1 = T0 + a1 * (h1 - h0)
                P1 = P0/((T1/T0)**(g/((a1/1000)*R)))
                rho1 = rho0/((T0/T1)**(-g/((a1/1000)*R) -1))
    
                Tprev = T0
                Pprev = P0
                rhoprev = rho0
    
            if i == 2:
                T2 = T1 + a2 * (h2 - h1)
                #ISOTHERMAL
                P2 = P1/(math.exp((g/(R*T2))*(h2*1000 - h1*1000)))
                rho2 = rho1*(math.exp((-g/(R*T2))*(h2*1000 - h1*1000)))
                
                Tprev = T1
                Pprev = P1
                rhoprev = rho1
    
            if i == 3:
                T3 = T2 + a3 * (h3 - h2)
                P3 = P2/((T3/T2)**( g/((a3/1000)*R)))
                rho3 = rho2/((T2/T3)**(-g/((a3/1000)*R) - 1))
                
                Tprev = T2
                Pprev = P2
                rhoprev = rho2
    
            if i == 4:
                T4 = T3 + a4 * (h4 - h3)
                P4 = P3/((T4/T3)**( g/((a4/1000)*R)))
                rho4 = rho3/((T3/T4)**(-g/((a4/1000)*R) - 1))
                
                Tprev = T3
                Pprev = P3
                rhoprev = rho3
                
            if i == 5:
                T5 = T4 + a5 * (h5 - h4)
                #ISOTHERMAL
                P5 = P4/(math.exp((g/(R*T5))*(h5*1000 - h4*1000)))
                rho5 = rho4*(math.exp((-g/(R*T5))*(h5*1000 - h4*1000)))
                
                Tprev = T4
                Pprev = P4
                rhoprev = rho4
    
            if i == 6:
                T6 = T5 + a6 * (h6 - h5)
                P6 = P5/((T6/T5)**( g/((a6/1000)*R)))
                rho6 = rho5/((T5/T6)**(-g/((a6/1000)*R) - 1))
                
                Tprev = T5
                Pprev = P5
                rhoprev = rho5
    
            if i == 7:
                T7 = T6 + a7 * (h7 - h6)
                P7 = P6/((T7/T6)**( g/((a7/1000)*R)))
                rho7 = rho6/((T6/T7)**(-g/((a7/1000)*R) - 1))
                
                Tprev = T6
                Pprev = P6
                rhoprev = rho6
                
            if i == zonelist[-1]:
                Tres = Tprev + ares * (X - hprev)
                
                if zonelist[-1] == 2 or zonelist[-1] == 5:
                
                    Pres = Pprev/(math.exp((g/(R*Tres))*(X*1000 - hprev*1000)))
                    rhores = rhoprev*(math.exp((-g/(R*Tres))*(X*1000 - hprev*1000)))
                
                else:
    
                    Pres = Pprev/((Tres/Tprev)**(g/((ares/1000)*R)))
                    rhores = rhoprev/((Tprev/Tres)**(- g/((ares/1000)*R) - 1))
        Tdifference = T0 - Tres
        Pper = (Pres / P0) * 100
        rhoper = (rhores / rho0) * 100

    elif check == False:
        X = 'NaN'
    elif X > 84 or X <0:
        X = -1
        
frame = simplegui.create_frame("Atmospheric pressure and air density calculator", 500, 300)
frame.set_canvas_background('Silver')
frame.add_input("Enter altitude in [km]:", inp, 100)
frame.set_draw_handler(draw)
frame.start()