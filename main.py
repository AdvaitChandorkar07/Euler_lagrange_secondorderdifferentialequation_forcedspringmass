import math
import matplotlib.pyplot as plt
import numpy as np
import csv

def func(t):
    '''This is a function that is used to calculate f(t)
    and it returns the value of f(t)
    for the time we enter'''
    #ret=math.cos(t)+math.sin(t)+1/2*math.pow(math.e,t)
    ret=math.cos(t)
    return ret

def z_dashcal(t,x,z):
    ret=(func(t)-k*x-c*z)/m
    return ret

def csv_filecreater(time,x):
    '''this function writes the value of the function
    at every interval in a csv file called csv which one can access
    '''
    with open('Data.csv', 'w') as csvfile:
        fieldnames = ['time', 'value']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        for i in range(len(time)):
            writer.writerow({'time': time[i], 'value': x[i]})


def graph_plotter(x,z,color,labe):
    '''Function used to plot the graph to
    display the solution of euler method '''
    plt.title('Differential Equation using euler method')
    #naming the x-axis and y-axis
    plt.xlabel('time')
    plt.ylabel('x(t)')
    #plt.plot(x, z)
    plt.plot(x, z, color=color, linewidth=1,  marker='o', markerfacecolor='blue', markersize=2,label=labe)

    #used to print all the calculated values
    #for i in range(int(partition)):
        #print(time[i], "        ", x[i])


#taking the input here
print("For the equation mx ̈ + cx ̇ + kx = f(t)")
m=float(input("value of m is ="))
c=float(input("value of c is ="))
k=float(input("value of k is ="))
delta=float(input("time step is ="))
ti=float(input("Initial time ="))
tf=float(input("Final time ="))
alpha=float(input("Enter alpha:-"))
#calculating all the extra values here
beta=1-alpha
partition=abs(ti-tf)/delta
no_of_terms=int((tf-ti)/delta+1)

i=ti
time=[]
while i<tf+delta:
    time.append((i))
    #time=[i for i in range(ti,tf+delta,delta)]
    i+=delta
#initializing the lists
x= [0 for i in range(len(time))]
z= [0 for i in range(len(time))]


#taking the initial values as entry
print("enter initial values ")
x[0]=float(input("for x initial"))
z[0]=float(input("for x' initial"))



def solver(alpha,beta,color,labe):
    '''this function where the equation is actually solved
    in a loop in iterative way
    and here we call the graph_plotter and csv_filecreater'''

    alpha,beta=beta,alpha
    # solving the main problem here
    for i in range(1,len(time)):
        A = np.array([[1, -1 * beta * delta], [beta * delta * k / m, 1 + c * beta * delta / m]])
        B = np.array([x[i - 1] + alpha * delta * z[i - 1],
                      z[i - 1] + alpha * delta * z_dashcal(time[i - 1], x[i - 1], z[i - 1]) + beta * func(
                          time[i]) * delta / m])
        [x[i], z[i]] = np.dot(np.linalg.inv(A), B)
    csv_filecreater(x, time)
    # calling the function to plot the data
    graph_plotter(time, x,color,labe)

#tp is the analytical solution which is used to plot
#tp=[(5/4*math.pow(math.e,i)- 1/4*math.pow(math.e,-i)-1/2*math.cos(i)) for i in time]
#plt.plot(time, tp, color='green', linewidth=1,  marker='o', markerfacecolor='blue', markersize=2,label='actual function')
# plotting for different values of alpha
solver(alpha, beta, 'blue', 'Euler Method')
#solver(0, 1, 'blue', 'alpha=1')
#solver(0.5, 0.5, 'black', 'alpha=0.5')
#solver(1, 0,'red','alpha=0')


plt.legend()
plt.show()










