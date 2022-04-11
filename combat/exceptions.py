class ConfoundingVariablesError(Exception):
    """Exception raised when confounding variables are detected.

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)
