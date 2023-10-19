import pickle
import random
import numpy as np
import concurrent.futures

# Definir el tamaño del sistema
tamaño_x = 16
tamaño_y = 16
tamaño = tamaño_x*tamaño_y

# Calcular la energía de un solo spin
def energySpin(matrix, x, y):
    rows, cols = matrix.shape
    spin_center = matrix[x, y]
    indices = [((x - 1) % rows, y), ((x + 1) % rows, y), (x, (y - 1) % cols), (x, (y + 1) % cols)]
    E = 0
    for i, j in indices:
        E += spin_center * matrix[i, j]
    return -E
                
# Calcular la magnetización de un estado (por spin)
def magnetization(matrix):
  sumatory = np.sum(matrix)
  N = tamaño
  M = np.abs((1/N)*sumatory)
  return M

# Caluclar la energía de un estado con matriz kernel y convoluciones para mayor eficiencia

def energiaEstado(matrix):
    def energy_convolution(matrix):
        kernel = np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]])
        conv_result = np.convolve(matrix.flatten(), kernel.flatten(), mode='same').reshape(matrix.shape)
        return -np.sum(conv_result * matrix)
    return energy_convolution(matrix)

# Algoritmo de metropolis

def algoritmo_metropolis(secuencia_inicial, temperatura, num_pasos,op):
    secuencia = [secuencia_inicial]
    estado_actual = secuencia_inicial
    energia = [energiaEstado(estado_actual)]
    magnetizacion = [magnetization(estado_actual)]
    for _ in range(num_pasos):
      for _ in range(1,tamaño):
        i, j = random.randint(0, estado_actual.shape[0] - 1), random.randint(0, estado_actual.shape[1] - 1)
        delta_energia = -2 * energySpin(estado_actual,i,j)
        # Distribución de boltzman
        if delta_energia < 0 or random.random() < np.exp(-delta_energia / temperatura): 
            estado_actual[i, j] *= -1
      secuencia.append(estado_actual)
      magnetizacion.append(magnetization(estado_actual))
      energia.append(energiaEstado(estado_actual))
    if op == 1:
      varE = np.var(energia[200:])
      varM = np.var(magnetizacion[200:])
      energia = np.mean(energia[200:])
      magnetizacion = np.mean(magnetizacion[200:])
    elif op==2:
      varE=np.var(energia)
      varM=np.var(magnetizacion)
    return secuencia,energia,magnetizacion,varE,varM

# Metodo para multiproccesing
def run_algoritmo_metropolis(configuracion_inicial,temperatura, num_pasos,op):
    return algoritmo_metropolis(configuracion_inicial, temperatura, num_pasos,op)

def testTemperature(configuracion_inicial,temperatura,op):
    num_iterations = 16
    num_pasos = 10000
    secuencias = []
    energias=[]
    magnetizaciones = []
    varEList = []
    varMList = []

    if __name__ == '__main__':
        with concurrent.futures.ProcessPoolExecutor() as executor:
            futures = [executor.submit(run_algoritmo_metropolis,configuracion_inicial, temperatura, num_pasos,op) for _ in range(num_iterations)]
            for future in concurrent.futures.as_completed(futures):
                secuencia_resultante,energia_resultante,magnetizacion_resultante,varE,varM = future.result()
                secuencias.append(secuencia_resultante)
                energias.append(energia_resultante)
                magnetizaciones.append(magnetizacion_resultante)
                varEList.append(varE)
                varMList.append(varM)
        return secuencias,energias,magnetizaciones,varEList,varMList

if __name__ == '__main__':
  op = 1 
  # Opción para correr una temperatura individual
  if op == 1:
    sizex = tamaño_x
    sizey = tamaño_y
    configuracion_inicial = np.random.choice([-1,1], size=(sizex, sizey))
    energias = []
    magnetizaciones = []
    varEnergias = []
    varMagnetizacion = []
    temperatures = np.linspace(1.8, 3, num=100)
    for i in range(len(temperatures)):
      secuencias,energia,magnetizacion,varE,varM = testTemperature(configuracion_inicial,temperatures[i],op)
      print("Temperature ",i,": Done!")
      with open("secuencia_estados"+str(i), "wb") as fp:
          pickle.dump(secuencias, fp)
      varEnergias.append(np.mean(varE))
      varMagnetizacion.append(np.mean(varM))
      energias.append(np.mean(energia))
      magnetizaciones.append(np.mean(magnetizacion))
    with open("energias", "wb") as fp:   
            pickle.dump(energias, fp)
    with open("magnetizaciones", "wb") as fp:   
            pickle.dump(magnetizaciones, fp)
    with open("varE", "wb") as fp:   
            pickle.dump(varEnergias, fp)
    with open("varM", "wb") as fp:   
            pickle.dump(varMagnetizacion, fp)
  elif op == 2:
    # Opción para correr una un rango de temperaturas
    sizex = tamaño_x
    sizey = tamaño_y
    configuracion_inicial = np.random.choice([-1,1], size=(sizex, sizey))
    temperature = 3
    secuencias,energia,magnetizacion,varEList,varMList = testTemperature(configuracion_inicial,temperature,op)
    for i in range(len(energia)):
      print(np.mean(energia[i]))
    with open("secuencia_estados3-0", "wb") as fp:   
      pickle.dump(secuencias, fp)
    with open("energias3-0", "wb") as fp:   
      pickle.dump(np.mean(energia,axis=0), fp)
    with open("magnetizaciones3-0", "wb") as fp:   
      pickle.dump(np.mean(magnetizacion,axis=0), fp)
    
