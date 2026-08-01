import numpy as np

class ConvergenceJudge:

    def __init__(self,
                 TargetState,
                 energies,
                 Values_ForceDisp):


        self.TargetState = TargetState
        self.energies = np.array(energies)
        self.NStates = self.energies.shape[1]
        self.TargetStateEnergies = np.array(self.energies[:,self.TargetState])
        self.max_force = np.array(Values_ForceDisp[0])
        self.rms_force = np.array(Values_ForceDisp[1])
        self.max_disp  = np.array(Values_ForceDisp[2])
        self.rms_disp  = np.array(Values_ForceDisp[3])
        self.threshold_maxforce  = float(Values_ForceDisp[4])
        self.threshold_rmsforce = float(Values_ForceDisp[5])
        self.threshold_maxdisp   = float(Values_ForceDisp[6])
        self.threshold_rmsdisp   = float(Values_ForceDisp[7])


    def StateCross(self):
        if self.TargetState > 0:

            States_X = []
            Ediff_withTarget = self.energies - self.TargetStateEnergies[:,None]
            for i in range(self.NStates):
                tmp_diff = Ediff_withTarget[:,i]
                idx = np.where(np.sign(tmp_diff[:-1]) != np.sign(tmp_diff[1:]))[0]
                if len(idx) != 0:
                    States_X.append(i)
    

            if States_X != []:
                
                return (
                    "transition",
                    0.1,
                    f"Crossing between {self.TargetState} and {','.join(map(str, States_X))}."
                )
            else:
                return 0
                    
        else:
            return 0


    def judge(self):

        diffs = (
            self.max_force[-1] -  self.threshold_maxforce,            
            self.rms_force[-1] -  self.threshold_rmsforce, 
            self.max_disp[-1]  -  self.threshold_maxdisp,   
            self.rms_disp[-1]  -  self.threshold_rmsdisp   
        )

        #
        #0. Conerged
        #
        if all(x < 0 for x in diffs):

            return (
                'converged', 
                 1.0,
                'Conerged!'
                )

        n = len(self.TargetStateEnergies)

        if n < 10:
            return (
                "unknown",
                0.5,
                "Too few optimization steps."
            )

        recent_E = self.TargetStateEnergies[-5:]
        recent_F = self.max_force[-5:]
        recent_D = self.max_disp[-5:]

        dE = np.diff(recent_E)

        abs_dE = np.abs(np.diff(recent_E))
        mean_dE = np.mean(abs_dE)

        force_ratio = recent_F[-1] / recent_F[0]
        disp_ratio = recent_D[-1] / recent_D[0]

        # Geometry optimizations do not necessarily improve monotonically.
        # Allow one temporary increase in the five-step window and regard the
        # displacement as decreasing only when the total decrease is meaningful.
        disp_decrease_fraction = np.mean(np.diff(recent_D) < 0)
        displacement_is_decreasing = (
            disp_ratio < 0.7 and disp_decrease_fraction >= 0.75
        )


        if self.NStates > 1:
           PossibleCIX =  self.StateCross()
           if PossibleCIX != 0:
                return PossibleCIX


        #
        # 1. Good progress
        #
        # Check this before stagnation: close to a minimum, both energy and
        # force can already be small while the optimizer is still reducing
        # the displacement.
        #
        if force_ratio < 0.5 and displacement_is_decreasing:

            return (
                "converging",
                0.8,
                "Force and displacement show sustained decrease."
            )

        #
        # 2. Stagnation
        #
        # A small energy change alone is not stagnation.  Only report it when
        # the displacement has failed to make meaningful recent progress.
        if (recent_F[-1] < 2e-4 and recent_D[-1] > 1e-3
                and not displacement_is_decreasing):

            if mean_dE < 1e-5:

                return (
                    "stagnation",
                    0.1,
                    "Force is small but displacement remains large."
                )

        #
        # 3. Diverse
        #
        if force_ratio > 2.0 and np.any(np.sign(dE[:-1]) != np.sign(dE[1:])) :

            return (
                "divergence",
                0.0,
                "Maximum force is increasing."
            )

        #
        # 4. Vibration
        #
        if np.std(recent_E) > 1e-3 and force_ratio > 0.8:

            return (
                "oscillation",
                0.3,
                "Energy oscillation detected."
            )

        #
        # 5. Almost converged
        #
        if recent_F[-1] < 5e-4 and recent_D[-1] < 5e-4:

            return (
                "near_convergence",
                0.95,
                "Near convergence."
            )

        return (
            "uncertain",
            0.5,
            "No clear trend."
        )
