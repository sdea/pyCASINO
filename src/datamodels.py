# Here we save the struct of the possible data models 
import numpy as np

class Pos:
    def __init__(self, rx, ry) -> None:
        self.rx = rx
        self.ry = ry

    def print_values(self):
        print(f'rx: {self.rx}')
        print(f'ry: {self.ry}')