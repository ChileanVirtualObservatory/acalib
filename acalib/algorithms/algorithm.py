class Algorithm(object):
    """
    Parent class of algorithms module. It defines the common
    interface that all algorithms should have.
    """
    def __init__(self, params=None):
        self.config = dict()
        if params is not None:
            for key,value in params.items():
                self.config[key] = value
        self.default_params()

    def default_params(self):
        pass

    def set_param(self, key, value):
        self.config[key] = value

    def get_param(self, key):
        if key in self.config:
            return self.config[key]
        else: return None

    def get_params(self):
        return self.config