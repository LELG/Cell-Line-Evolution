import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.path import Path
from matplotlib.spines import Spine
from matplotlib.projections.polar import PolarAxes
from matplotlib.projections import register_projection

def _radar_factory(num_vars):
    theta = 2*np.pi * np.linspace(0, 1-1./num_vars, num_vars)
    theta += np.pi/2

    def unit_poly_verts(theta):
        x0, y0, r = [0.5] * 3
        verts = [(r*np.cos(t) + x0, r*np.sin(t) + y0) for t in theta]
        return verts

    class RadarAxes(PolarAxes):
        name = 'radar'
        RESOLUTION = 1

        def fill(self, *args, **kwargs):
            closed = kwargs.pop('closed', True)
            return super(RadarAxes, self).fill(closed=closed, *args, **kwargs)

        def plot(self, *args, **kwargs):
            lines = super(RadarAxes, self).plot(*args, **kwargs)
            for line in lines:
                self._close_line(line)

        def _close_line(self, line):
            x, y = line.get_data()
            # FIXME: markers at x[0], y[0] get doubled-up
            if x[0] != x[-1]:
                x = np.concatenate((x, [x[0]]))
                y = np.concatenate((y, [y[0]]))
                line.set_data(x, y)

        def set_varlabels(self, labels):
            self.set_thetagrids(theta * 180/np.pi, labels)

        def _gen_axes_patch(self):
            verts = unit_poly_verts(theta)
            return plt.Polygon(verts, closed=True, edgecolor='k')

        def _gen_axes_spines(self):
            spine_type = 'circle'
            verts = unit_poly_verts(theta)
            verts.append(verts[0])
            path = Path(verts)
            spine = Spine(self, spine_type, path)
            spine.set_transform(self.transAxes)
            return {'polar': spine}

    register_projection(RadarAxes)
    return theta

def radar_graph(labels = [], values = [], optimum=[]):
    N = len(labels) 
    theta = _radar_factory(N)
    max_val = max(max(optimum), max(values))
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection='radar')
    ax.plot(theta, values, color='k')
    ax.plot(theta, [1,2,3,4,3], color = 'k')
    ax.plot(theta, optimum, color='r')
    ax.set_varlabels(labels)
    plt.show()
    #plt.savefig("radar.png", dpi=100)

def generate_Radar_plot(labels, BestFit, WorstFit):
	N = len(labels) 
	theta = _radar_factory(N)
	max_val = max(max(BestFit), max(WorstFit))
	fig = plt.figure()
	ax = fig.add_subplot(1, 1, 1, projection='radar')
	for value in BestFit:
		ax.plot(theta, value, color='y')
	for value in WorstFit:
		ax.plot(theta, value, color='b')
	ax.set_varlabels(labels)
	plt.show()

def generate_Single_Radar_plot(labels, Data):
	N = len(labels)
	theta = _radar_factory(N)
	max_val = max(Data)
	fig = plt.figure()
	fig.suptitle('Bestfit Geneation', fontsize=14, fontweight='bold')
	ax = fig.add_subplot(1, 1, 1, projection='radar')
	in_ = True
	for value in Data:
		if in_:
			in_ = False
			ax.plot(theta, value, color='b', label = 'Best Fit')
		else:
			ax.plot(theta, value, color = 'b')
	ax.set_varlabels(labels)# put the labe up
	ax.legend()
	plt.show()

def generate_Single_BestFitB4(bestFit, generation, save_file_path):
	labels = ['P_Ben', 'P_Del', 'P_Inc', 'P_Dec', 'Fitness']
	N = len(labels)
	theta = _radar_factory(N)
	max_val = max(bestFit)
	fig = plt.figure()
	title = 'Best Memebers In Geneation '+ generation
	fig.suptitle(title, fontsize=14, fontweight='bold')
	ax = fig.add_subplot(1, 1, 1, projection='radar')
	in_ = True
	for element in bestFit:
		if in_:
			in_ = False
			ax.plot(theta, element, color='b', label = 'Best Fit')
		else:
			ax.plot(theta, element, color = 'b')
	ax.set_varlabels(labels)
	ax.legend()# put the label in ax.plot
	#plt.show()
	plt.savefig(save_file_path)

def generate_Best_and_Worst(BestFit, WorstFit, generation, save_file_path):
	labels = ['P_Ben', 'P_Del', 'P_Inc', 'P_Dec', 'Fitness']
	N = len(labels)
	theta = _radar_factory(N)
	max_val = max(max(BestFit), max(WorstFit))
	fig = plt.figure()
	title = 'Best and Worst Memebers In Geneation '+ generation
	fig.suptitle(title, fontsize=14, fontweight='bold')
	ax = fig.add_subplot(1, 1, 1, projection='radar')
	_in = True
	for best, worst in zip(BestFit, WorstFit):
		if _in:
			_in = False
			ax.plot(theta, best, color = 'b', label = 'Best Fit')
			ax.plot(theta, worst, color = 'y', label = 'Worst Fit')
		else:
			ax.plot(theta, best, color = 'b')
			ax.plot(theta, worst, color = 'y')
	ax.set_varlabels(labels)
	ax.legend()# put the label in ax.plot
	#plt.show()
	plt.savefig(save_file_path)






if __name__ == '__main__':
	#tester framework
 	labels = ['P_Ben', 'P_Del', 'P_Inc', 'P_Dec', 'Fitness']
 	values =  [[1, 1, 2, 7, 4],
 	           [1, 2, 3, 4, 5]]
 	optimum = [[5, 3, 2, 4, 5],
 	           [1, 2, 3, 4, 5]]
 	generate_Radar_plot(labels, values, optimum)
 	generate_Single_Radar_plot(labels, values)