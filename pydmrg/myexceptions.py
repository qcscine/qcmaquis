# new exception classes for potential catching
# at the moment, they just print text via the Exception class

class ConstantError(Exception):
    pass

class MPSTensorError(Exception):
    pass

class MPOTensorError(Exception):
    pass

class SingularValueTruncationError(Exception):
    pass

class DecompositionError(Exception):
    pass

class DimensionMismatchError(Exception):
    pass

class MatrixProductStateDimensionError(Exception):
    pass

class MatrixProductOperatorDimensionError(Exception):
    pass

class NormalizationError(Exception):
    pass

class CanonicalFormMismatchError(Exception):
    pass

class StateInitializationError(Exception):
    pass

