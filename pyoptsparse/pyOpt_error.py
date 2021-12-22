"""
pyOptSparse_error

Holds a simple error handling class for pyOptSparse
"""


class Error(Exception):
    """
    Format the error message in a box to make it clear this
    was a expliclty raised exception.
    """

    def __init__(self, message):
        msg = "\n+" + "-" * 78 + "+" + "\n" + "| pyOptSparse Error: "
        i = 21
        for word in message.split():
            if len(word) + i + 1 > 78:  # Finish line and start new one
                msg += " " * (79 - i) + "|\n| " + word + " "
                i = 2 + len(word) + 1
            else:
                msg += word + " "
                i += len(word) + 1
        msg += " " * (79 - i) + "|\n" + "+" + "-" * 78 + "+" + "\n"
        print(msg)
        self.message = message
        super().__init__(message)

    def __str__(self):
        return self.message


class pyOptSparseWarning:
    """
    Format a warning message
    """

    def __init__(self, message):
        msg = "\n+" + "-" * 78 + "+" + "\n" + "| pyOptSparse Warning: "
        i = 23
        for word in message.split():
            if len(word) + i + 1 > 78:  # Finish line and start new one
                msg += " " * (79 - i) + "|\n| " + word + " "
                i = 1 + len(word) + 1
            else:
                msg += word + " "
                i += len(word) + 1
        msg += " " * (78 - i) + "|\n" + "+" + "-" * 78 + "+" + "\n"
        print(msg)
