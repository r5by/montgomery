from abc import ABC, abstractmethod


class GFType(ABC):

    def __init__(self, order):
        self.order = order
    @abstractmethod
    def __call__(self, value) -> 'GFElementType':
        """ Converts value into this field"""
        pass


class GFElementType(ABC):
    @abstractmethod
    def __add__(self, other: 'GFElementType') -> 'GFElementType':
        pass

    @abstractmethod
    def __sub__(self, other: 'GFElementType') -> 'GFElementType':
        pass

    @abstractmethod
    def __mul__(self, other: 'GFElementType') -> 'GFElementType':
        pass

    @abstractmethod
    def __truediv__(self, other: 'GFElementType') -> 'GFElementType':
        pass

    @abstractmethod
    def __rtruediv__(self, other: int) -> 'GFElementType':
        pass

    @abstractmethod
    def __neg__(self) -> 'GFElementType':
        pass

    @abstractmethod
    def __pow__(self, exponent: int) -> 'GFElementType':
        pass

    @abstractmethod
    def __int__(self) -> int:
        pass

    @abstractmethod
    def __eq__(self, other: 'GFElementType') -> bool:
        pass
