# 绘图部分
import matplotlib.pyplot as plt


def CrackLength_Cycles(sequence, n_MTS, a_MTS, n_Manual=[], a_Manual=[], show=0, save=0):
    # Plot: Crack Length - Cycles by Manual and MTS
    plt.figure(num=1, figsize=(7, 5))
    plt.scatter(n_MTS, a_MTS, s=1, label='$MTS$', color='red', lw=1)
    if n_Manual and a_Manual:
        plt.scatter(n_Manual, a_Manual, s=1, label='$Manual$', color='blue', lw=1)
    plt.xlabel("Cycles")
    plt.ylabel("Crack Length (mm)")
    plt.title('Crack Length - Cycles '+sequence+'(MTS compare with Manual)')
    plt.legend()
    plt.grid()
    if save:
        plt.savefig('cracklength_cycels'+sequence+'.png', dpi=320)
    if show:
        plt.show()


def CrackLengthError_Cycles(sequence, a_Manual, a_MTS, n, show=0, save=0):
    # Plot: Crack Length - Cycles Error by MTS and Manual
    plt.figure(num=2, figsize=(7, 5))
    plt.scatter(n, a_Manual - a_MTS, s=1, label='$Manual - MTS$', color='black', lw=1)
    plt.xlabel("Cycles")
    plt.ylabel("Crack Length Error(mm)")
    plt.title('Crack Length - Cycles Error'+sequence+'(MTS compare with Manual)')
    plt.legend()
    plt.grid()
    if save:
        plt.savefig('cracklength_cycels_error'+sequence+'.png', dpi=320)
    if show:
        plt.show()


def FCGR_DeltaK_MTS(sequence, dk, dadn, dadn_paris=[], dadn_walker=[], dadn_reference=[], show=0, save=0):
    # Plot: da/dN - dk plot(MTS Result)
    plt.figure(num=3, figsize=(7,5))
    plt.scatter(dk, dadn, s=1, label='$Experiment$', color='red', lw=1)
    if len(dadn_paris):
        plt.plot(dk, dadn_paris, label='$Fitting By Paris$', color='blue', linewidth=2)
    if len(dadn_walker):
        plt.plot(dk, dadn_walker, label="$Walker's Model$", color='purple', linewidth=2)
    if len(dadn_reference):
        plt.plot(dk, dadn_reference, label="$Reference Paris$", color='black', linewidth=2)
    plt.xlabel("DeltaK Applied (MPa*m^0.5)")
    plt.xscale('log')
    plt.ylabel("da/dN (mm/cycle)")
    plt.yscale('log')
    plt.title('da/dN - dK '+sequence+'(MTS Result)')
    plt.legend()
    plt.grid()
    if save:
        plt.savefig('dadn_dk_MTS'+sequence+'.png', dpi=320)
    if show:
        plt.show()


def FCGR_DeltaK_Comparation(sequence, dk_MTS, dadn_MTS, dk_Manual, dadn_manual, dadn_paris_MTS = [], dadn_paris_Manual = [], show=0, save=1):
    # Plot: da/dN - dk plot(MTS compare with Manual)
    plt.figure(num=4,figsize=(7,5))
    plt.scatter(dk_MTS, dadn_MTS, s=1, label='$MTS Results$', color='red', lw=1)
    plt.scatter(dk_Manual, dadn_manual, s=1, label='$Manual Results$', color='blue', lw=1)
    if len(dadn_paris_MTS):
        plt.plot(dk_MTS, dadn_paris_MTS, label='$MTS Results Fitting$', color='red', linewidth=2)
    if len(dadn_paris_Manual):
        plt.plot(dk_Manual, dadn_paris_Manual, label='$Manual Results Fitting$', color='blue', linewidth=2)
    plt.xlabel("DeltaK Applied (MPa*m^0.5)")
    plt.xscale('log')
    plt.ylim(min(dadn_MTS)*0.9, max(dadn_MTS)*1.1)
    plt.xlim(min(dk_MTS)*0.9, max(dk_MTS)*1.1)
    plt.ylabel("da/dN (mm/cycle)")
    plt.yscale('log')
    plt.title('da/dN - dK '+sequence+'(MTS compare with Manual)')
    plt.legend()
    plt.grid()
    if save:
        plt.savefig('dadn_dk_MTS_and_Manual'+sequence+'.png', dpi=320)
    if show:
        plt.show()
