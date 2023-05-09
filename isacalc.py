# Python ISA calculator by Florian Voetter 5326168
import math as math
import re

#constants and conversion factors
ft = 0.3048
FL = 100.0 * ft
g_0 = 9.80665
R = 287.0
T_SL = 288.15
p_SL = 101325.0


#Order units as a dictionary
units = {
	'1':{'name':'meters', 'conv':1.0},
	'2':{'name':'feet', 'conv':ft},
	'3':{'name':'FL', 'conv':FL},
	'4':{'name':'press'},
	'5':{'name':'dens'}
	}


#Make an array with the properties of each layer.

#First define the class
class layer:
	def __init__(self, name, h_b, h_e, a):
		self.name = name
		self.h_b = h_b
		self.h_e = h_e
		self.a = a
		
#Then enter the data and order them insde an array. This makes it easy to say implement another layer
layers = [
	layer('Troposphere', 0, 11000, -0.0065),
	layer('Tropopause', 11000, 20000, 0),
	layer('Stratosphere 1', 20000, 32000, 0.0010),
	layer('Stratosphere 2', 32000, 47000, 0.0028),
	layer('Stratopause', 47000, 51000, 0),
	layer('Mesophere 1', 51000, 71000, -0.0028),
	layer('Mesophere 2', 71000, 89000, -0.0020)
	]

#Define functions
def rho_at(p, T):
	return p/(R*T)
	
	
def T_at(T_0, dh, a):
	return T_0+a*dh

def p_grad(p_0, T_0, T_1, a):
	return p_0*pow((T_1/T_0),-(g_0/(a*R)))

def p_isot(p_0, T, dh):
	return  p_0*math.exp((-g_0*dh)/(R*T))

#Inverse version of p_grad
def dh_grad(p_1, p_0, T_0, a):
	return (T_0/a)*(pow((p_1/p_0),((-a*R)/g_0))-1)
	
#Inverse version of p_isot
def dh_isot(p_1, p_0, T):
	return ((-R*T)/g_0)*math.log(p_1/p_0)

def isa(h, p_sl=p_SL, T_sl=T_SL):
	h_cur=0.0
	p_cur=p_sl
	T_base=T_sl
	T_cur=T_sl
	#The algorithm itself consist of a for loop which loops through each layer. 
	for layer_cur in layers:
		#Then the height delta for the layer is calculated. If the input height is higher than the end height of the layer, the height delta is just the height of the current layer.
		#If the input height is between the end and start of a layer, then the height to calculate over is the difference between the input height and the base height of the current layer.
		#If the input height is lower than the base height of layer, the height delta is zero. Thus the pressure won't be changed.
		dh=max(0,min(h-layer_cur.h_b,layer_cur.h_e-layer_cur.h_b))
		if  layer_cur.a != 0: 
			#Calc pressure drop if layer is a gradient layer. First the temperature at the top is calculated, the the pressure drop.
			T_cur = T_at(T_base, dh, layer_cur.a)
			p_cur = p_grad(p_cur, T_base, T_cur, layer_cur.a)
			#The temperature at the top of the current layer will be the temperature at the base of the next layer.
			T_base = T_cur
		else:
			#Calc pressure drop if layer is a isotherm layer. 
			p_cur = p_isot(p_cur, T_cur, dh)
	return {
		'press':p_cur,
		'temp':T_cur,
		'dens':rho_at(p_cur, T_cur)
		}

def isa_inv_p(p, p_sl, T_sl):
	
	dh=0
	h=0
	#Find in which layer the value of p is
	for layer_cur in layers:
		if p >= isa(layer_cur.h_e, p_sl, T_sl)['press']:
			isa_b=isa(layer_cur.h_b, p_sl, T_sl)

			if  layer_cur.a != 0: 
				#Gradient layer
				dh=dh_grad(p, isa_b['press'], isa_b['temp'], layer_cur.a)
			else:
				#Isotherm layer
				dh=dh_isot(p, isa_b['press'], isa_b['temp'])
			
			h=layer_cur.h_b+dh
			break
		
	return h

def isa_inv_rho(rho, p_sl, T_sl):
	
	dh=0
	h=0
	#Find in which layer the value of rho is
	for layer_cur in layers:
		if rho >= isa(layer_cur.h_e, p_sl, T_sl)['rho']:
			isa_b=isa(layer_cur.h_b, p_sl, T_sl)

			if  layer_cur.a != 0: 
				#Gradient layer
				dh=dh_grad(p, isa_b['press'], isa_b['temp'], layer_cur.a)
			else:
				#Isotherm layer
				dh=dh_isot(p, isa_b['press'], isa_b['temp'])
			
			h=layer_cur.h_b+dh
			break
		
	return h	


	


if __name__ == '__main__':
	print(
		'**** ISA calculator ****\n'
		'1. Calculate ISA for altitude in meters\n'
		'2. Calculate ISA for altitude in feet\n' 
		'3. Calculate ISA for altitude in FL\n' 
		'4. Calculate height from given pressure \n'
		'5. Calculate height from given density \n'
		)



	#Allow the user to select their preffered unit. If the user types in 7 or something the program will just default to meters.
	selection = input('Enter your choice:')
	if selection.isdigit() == False or int(selection) < 1 or int(selection) > 5:
		selection='1'
		print('Unclear selection. Defaulting to selection 1.')

	print('Do you want to define a different sea level temperature than 288.15K?\n')
	def_slt=str(input('[y/n]:'))
	if def_slt=='y':
		T_SL=float(input('Enter new sea level temp[K]:'))

	if int(selection) <= 3:
		

		#Here the altitude is actually retrieved by calling get_h(), and the ValueError is wrapped. Note that ValueError will also be called if the user types a string instead of a number.
		h_in=0
		error = True
		while error: 
			try:
				h_in = float(input('Enter altitude[{unit}]:'.format(unit=units[selection]['name'])))
				print(h_in)
				if h_in < 0.0 or h_in > layers[-1].h_e/units[selection]['conv']:
					range_error = ValueError()
					raise range_error
			except ValueError:
				print('please input a valid number between 0 and {maxalt}'.format(maxalt=round(layers[-1].h_e/units[selection]['conv'])))
			else:
				error = False
				
		#Convert unit to meters
		h_in = h_in * units[selection]['conv']
		print(h_in)

		#Call the ISA algorithm to calculate values
		isa_final = isa(h_in, p_SL, T_SL)

		print('**** final values ****')	
		print('Temperature: {T_K} K, {T_C} \'C'.format(T_K = round(isa_final['temp']*100)/100, T_C = round((isa_final['temp']-273.15)*100)/100))
		print('Pressure: {press} Pa ({press_rat} %SL)'.format(press=round(isa_final['press']*100)/100, press_rat=round((isa_final['press']/p_SL)*100)))
		print('Density: {dens} kg/m^3 ({dens_rat} %SL)'.format(dens=round(isa_final['dens']*100000)/100000, dens_rat=round((isa_final['dens']/rho_at(p_SL, T_SL))*100)))


	else:
		error=True
		while error: 
			try:
				p_in=float(input('Enter pressure[Pa]:'))
				if p_in <= 0.0 or p_in > p_SL:
					range_error = ValueError()
					raise range_error
			except ValueError:
				print('please input a valid number between {minp} and {maxp}'.format(minp=math.ceil(10*isa(layers[-1].h_e, p_SL, T_SL)['press'])/10, maxp=p_SL))
			else:
				error = False
		
		
		h_final=isa_inv_p(p_in, p_SL, T_SL)
		
		print('Height: {h_m} m, {h_feet} ft, {h_FL} FL'.format(h_m=round(h_final), h_feet=round(h_final/ft), h_FL=round(h_final/FL)))
		
	dummy = input('Press the any key to end the ISA calculator.')	




	
	







