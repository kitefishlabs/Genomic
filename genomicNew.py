def Individual():
	
	def __init__(idnum, params):
	
	def analyze(idnum)
		returns np.ndarray
	def render(idnum, corpus, )
		returns sc.Synth



def Population():

	def __init(size):
	
	

def GeneSpec():

	def __init__(self, arg=None, **gspec_params):
		self._initialize(analysis_params)
	
	
	def _initialize(self, gspec_params=None):
		self._check_spec(gspec_params)
	
	def _check_spec(self, gspec_params=None)
	
		self.gspec_params = gspec_params if gspec_params is not None else self.gspec_params
		dgsp = self.default_gspec_params()
		for k in dgsp.keys():
			self.gspec_params[k] = self.gspec_params.get(k, dgsp[k])
		return self.gspec_params
			
	@staticmethod
	def default_gspec_params():
		gspec_params = {
			'clip' : 1.0,
			'amp0' : 1.0
		}
		return gspec_params
	
	def 