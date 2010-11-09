import matplotlib.pyplot as plt

def plotevals(data, folder):
    plt.figure()
    for d in data:
        plt.plot(d[:]['v'], d[:]['lambda1'], 'k,')
        plt.plot(d[:]['v'], d[:]['lambda2'], 'k,')
        plt.plot(d[:]['v'], d[:]['lambda3'], 'k,')
        plt.plot(d[:]['v'], d[:]['lambda4'], 'k,')
    # Presumably, all eigenvalues are calculated over same speed range
    plt.axis([min(d[:]['v']), max(d[:]['v']), -10, 10])
    plt.xlabel('Speed [m/s]')
    plt.grid(True)
    plt.savefig(folder + 'evals.pdf')

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

def plotcontact(d, folder):
    plt.figure()
    plt.plot(d[:]['q6'], d[:]['q7'], 'r-', label='Rear contact')
    plt.plot(d[:]['fnx'], d[:]['fny'], 'g-', label='Front contact')
    plt.xlabel('x [meters]')
    plt.ylabel('y [meters]')
    plt.axis('tight')
    plt.title('Contact point locations in x-y plane')
    plt.legend(loc=0)
    plt.savefig(folder + 'xy.pdf')

def plotenergy(d, folder):
    plt.figure()
    plt.plot(d[:]['t'], d[:]['ke'], 'r-', label='KE')
    plt.plot(d[:]['t'], d[:]['pe'], 'g-', label='PE')
    plt.plot(d[:]['t'], d[:]['pe'] + d[:]['ke'], 'b-', label='TE')
    plt.xlabel('seconds')
    plt.ylabel('Joules')
    plt.title('Energy')
    plt.legend(loc=0)
    plt.savefig(folder + 'energy.pdf')

def plotconstraintforces(d, folder):
    f2, (f2a1, f2a2, f2a3) = plt.subplots(3, sharex=True, sharey=False)
    f2a1.plot(d[:]['t'], d[:]['Rx'], 'r-', label='$R_x$')
    f2a1.plot(d[:]['t'], d[:]['Fx'], 'g-', label='$F_x$')
    f2a1.set_title('Wheel contact forces')
    f2a1.legend(loc=0)
    f2a1.set_yticks(f2a1.get_yticks()[1:])
    f2a2.plot(d[:]['t'], d[:]['Ry'], 'r-', label='$R_y$')
    f2a2.plot(d[:]['t'], d[:]['Fy'], 'g-', label='$F_y$')
    f2a2.legend(loc=0)
    f2a2.set_yticks(f2a2.get_yticks()[1:-1])
    f2a3.plot(d[:]['t'], d[:]['Rz'], 'r-', label='$R_z$')
    f2a3.plot(d[:]['t'], d[:]['Fz'], 'g-', label='$F_z$')
    f2a3.legend(loc=0)
    f2a3.set_yticks(f2a3.get_yticks()[:-1])
    f2a3.axes.set_xlabel('seconds')
    # Fine-tune figure; make subplots close to each other and hide x ticks for
    # all but bottom plot.
    f2.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in f2.axes[:-1]], visible=False)
    plt.savefig(folder + 'constraintforces.pdf')

def plotconstraints(d, folder):
    f3, (f3a1, f3a2, f3a3, f3a4, f3a5) = plt.subplots(5, sharex=True, sharey=False)
    f3a1.plot(d[:]['t'], d[:]['fnz'], label='Holonomic constraint')
    f3a1.set_title('Numerical satisfaction of constraints')
    f3a1.legend(loc=0)
    f3a1.set_yticks(f3a1.get_yticks()[1:])
    f3a2.plot(d[:]['t'], d[:]['nh1'], label='Non-holonomic constraint 1')
    f3a2.legend(loc=0)
    f3a2.set_yticks(f3a2.get_yticks()[1:-1])
    f3a3.plot(d[:]['t'], d[:]['nh2'], label='Non-holonomic constraint 2')
    f3a3.legend(loc=0)
    f3a3.set_yticks(f3a3.get_yticks()[1:-1])
    f3a4.plot(d[:]['t'], d[:]['nh3'], label='Non-holonomic constraint 3')
    f3a4.legend(loc=0)
    f3a4.set_yticks(f3a4.get_yticks()[1:-1])
    f3a5.plot(d[:]['t'], d[:]['ke'] + d[:]['pe'] -
                           (d[0]['ke'] + d[0]['pe']), label='$\Delta(KE+PE)$')
    f3a5.legend(loc=0)
    f3a5.set_yticks(f3a5.get_yticks()[:-1])
    f3.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in f3.axes[:-1]], visible=False)
    f3a5.axes.set_xlabel('seconds')
    plt.savefig(folder + 'constraints.pdf')

def timeseriesplots(plot_dict, sim_data, folder):
    if plot_dict['orientation']:
        plotorientation(sim_data, folder)
    if plot_dict['contact']:
        plotcontact(sim_data, folder)
    if plot_dict['energy']:
        plotenergy(sim_data, folder)
    if plot_dict['constraintforces']:
        plotconstraintforces(sim_data, folder)
    if plot_dict['constraints']:
        plotconstraints(sim_data, folder)

def plotcontroller(plot_dict,
                   sim_data=None,
                   eval_data=None,
                   steady_turning_data=None,
                   folder=None):

    if folder == None:
        folder = ""
    if sim_data != None:
        timeseriesplots(plot_dict, sim_data, folder)
    if eval_data != None and plot_dict['evals']:
        plotevals(eval_data, folder)
