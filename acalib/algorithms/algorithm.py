class Algorithm(object):
    """
    Parent class of algorithms module. It defines the common
    interface that all algorithms should have.
    """
    def __init__(self, params=None):
        """
        Load default params if None given.

        Parameters
        ----------            
            params : dict (default = None)
                Dictionary with algorithm's parameters
        """
        self.config = dict()
        if params is not None:
            for key,value in params.items():
                self.config[key] = value
        self.default_params()

    def default_params(self):
        """
            Set default params to algorithm's cofiguration. Every algorithm has to implement their own default parameters.
        """
        pass

    def set_param(self, key, value):
        """
            Set the value of a parameter given a key.
            
            Parameters
            ----------     
            key : string
                algorithm internal config dictionary key.
            value :  
                value to set.
        """
        self.config[key] = value

    def get_param(self, key):
        """
            Get the value of a parameter given a key.
            
            Parameters
            ----------     
            key : string
                algorithm internal config dictionary key.

            Returns
            -------
            value for the given key.            
        """
        if key in self.config:
            return self.config[key]
        else:
            return None

    def get_params(self):
        """
            Get the internal configuration dictionary.

            Returns
            -------
            dict with the current algorithm's configuration.
        """
        return self.config
