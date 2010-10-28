import matplotlib.pyplot as plt

def plotevals(d, folder):
    plt.figure()
    plt.plot(d[:]['v'], d[:]['lambda1'], 'k,')
    plt.plot(d[:]['v'], d[:]['lambda2'], 'k,')
    plt.plot(d[:]['v'], d[:]['lambda3'], 'k,')
    plt.plot(d[:]['v'], d[:]['lambda4'], 'k,')
    plt.axis([min(d[:]['v']), max(d[:]['v']), -10, 10])
    #plt.savefig(folder + 'evals.pdf')

def plotfunctions(plot_dict, data, folder):
    if plot_dict['evals']:
        plotevals(data, folder)

def plotorientation(d, folder):
    f1, (f1a1, f1a2, f1a3) = plt.subplots(3, sharex=True, sharey=False)
    f1a1.plot(d[:]['t'], d[:]['q0'], label='Frame Yaw')
    f1a1.plot(d[:]['t'], d[:]['q1'], label='Frame Lean')
    f1a1.plot(d[:]['t'], d[:]['q2'], label='Frame Pitch')
    f1a1.legend(loc=0)
    f1a1.set_title('Bicycle orientation angles')
    f1a1.set_yticks(f1a1.get_yticks()[1:])
    f1a2.plot(d[:]['t'], d[:]['q3'], label='Steer')
    f1a2.legend(loc=0)
    f1a2.set_yticks(f1a2.get_yticks()[1:-1])
    f1a3.plot(d[:]['t'], d[:]['fa_yaw'], label='Fork Yaw')
    f1a3.plot(d[:]['t'], d[:]['fa_lean'], label='Fork Lean')
    f1a3.plot(d[:]['t'], d[:]['fa_pitch'], label='Fork Pitch')
    f1a3.legend(loc=0)
    f1a3.set_yticks(f1a3.get_yticks()[:-1])
    f1.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in f1.axes[:-1]], visible=False)
    f1a3.axes.set_xlabel('seconds')
    plt.savefig(folder + 'orientation.pdf')
