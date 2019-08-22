import matplotlib.pyplot as plt


def plot_f(errors, values):
    plt.plot(values, errors)
    plt.ylabel("Error")
    plt.xlabel("Iteraciones")
    plt.show()


def plot(errors, title):
    plt.plot(errors)
    plt.title(title)
    plt.ylabel("Error")
    plt.xlabel("Iteraciones")
    plt.show()