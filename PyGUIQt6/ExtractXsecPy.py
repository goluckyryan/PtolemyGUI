#!/usr/bin/python3

import numpy as np

def extract_xsec(read_file, index_for_elastic=1):
    # index_for_elastic = 1 ; for Ratio
    # index_for_elastic = 2 ; for Total
    # index_for_elastic = 3 ; for Rutherford

    # --- open file
    try:
        with open(read_file, 'r') as file_in:
            # ---- variable for Xsec
            title = []
            data_matrix = []
            angle = []
            Ex = []
            total_xsec = []
            reaction = []
            angle_min = 0
            angle_max = 0
            angle_step = -1

            # ================================== read file
            line_num = 0
            data_xsec = []
            start_extract = False
            angle_filled = False
            num_cal = 0
            reaction_flag = 0
            pre_find_for_elastic = False
            print("======================================================")
            
            for line in file_in:
                line_num += 1
                line = line.strip()

                # ---- find Title
                find_str = "$============================================"
                pos = line.find(find_str)
                if pos != -1:
                    title.append(line[pos + len(find_str) + 1:])
                    pos = title[-1].find("(")
                    ex = title[-1][3:pos]
                    Ex.append(float(ex))
                    data_xsec.clear()
                    num_cal += 1
                    continue

                # ---- find Reaction
                find_str = "0INPUT... CHANNEL"
                pos = line.find(find_str)
                if pos != -1:
                    reaction.append(line[pos + len(find_str) + 1:])
                    reaction_flag = 1  # 1 for (d,d) or (p,p)
                    continue

                find_str = "REACTION:"
                pos = line.find(find_str)
                if pos != -1:
                    reaction.append(line[pos + len(find_str) + 1:])
                    reaction_flag = 2  # 2 for (d,p) or (p,d)
                    continue

                # ----- pre find for Elastic scattering
                find_str = "0        OPTICAL MODEL SCATTERING FOR THE OUTGOING CHANNEL"
                pos = line.find(find_str)
                if pos != -1:
                    pre_find_for_elastic = True
                    continue

                # ----- find angle setting when not known
                if angle_step == -1:
                    find_str = "anglemin="
                    pos = line.find(find_str)
                    if pos != -1:
                        pos1 = line.find("anglemax=")
                        angle_min = float(line[pos + len(find_str):pos1].strip())
                        pos = pos1
                        pos1 = line.find("anglestep=")
                        angle_max = float(line[pos + len(find_str):pos1].strip())
                        angle_step = float(line[pos1 + len(find_str) + 1:].strip())
                        continue

                # ----- check if start extracting Xsec or not 
                find_str = "dumpdumpdump"
                if reaction_flag == 1 and pre_find_for_elastic:
                    find_str = "C.M.    LAB     RUTHERFORD"
                if reaction_flag == 2:
                    find_str = "0  C.M.  REACTION     REACTION   LOW L  HIGH L   % FROM"
                pos = line.find(find_str)
                if pos != -1:
                    start_extract = True
                    continue

                # ----- start extracting Xsec
                if start_extract:
                    if len(line) < 20:
                        continue

                    num_angle = num_xsec = None
                    if reaction_flag == 1:  # Elastics (d,d) or (p,p)
                        num_angle = float(line[0:7].strip())
                        if index_for_elastic == 1:
                            num_xsec = float(line[15:30].strip())
                        elif index_for_elastic == 2:
                            if len(line) > 60:
                                num_xsec = float(line[30:43].strip())
                            else:
                                num_xsec = -404
                        else:
                            if len(line) > 60:
                                num_xsec = float(line[57:71].strip())
                            else:
                                num_xsec = -404

                    if reaction_flag == 2:  # InElastics (d,p) or (p,d)
                        if is_float(line[0:7].strip()):
                            num_angle = float(line[0:7].strip())
                            num_xsec = float(line[7:26].strip())
                        else:
                            num_angle = -1.0
                            num_xsec = -1.0

                    if num_angle != num_xsec and num_xsec > 0.:
                        if not angle_filled:
                            angle.append(num_angle)
                        data_xsec.append(num_xsec)

                # ------ find total Xsec, if found stop extraction
                find_str = "dumpdumpdump"
                if reaction_flag == 1 and pre_find_for_elastic:
                    find_str = "0TOTAL REACTION CROSS SECTION ="
                if reaction_flag == 2:
                    find_str = "0TOTAL:"
                pos = line.find(find_str)

                if pos != -1:
                    total_xsec.append(float(line[pos + len(find_str):].strip()))
                    print(f"{num_cal:2d} | {title[-1]} | total Xsec(4pi): {total_xsec[-1]} mb")
                    data_matrix.append(data_xsec[:])
                    angle_filled = True
                    start_extract = False
                    reaction_flag = 0
                    pre_find_for_elastic = False
                    continue

        # ================================== summary
        print("====================== Total number of lines read : ", line_num)
        print(f"Angle : {angle_min:5.2f}, {angle_max:5.2f} | step : {angle_step:5.2f}")
        print("Number of Angles read : ", len(angle))
        print("Number of data read : ", len(data_xsec))
        print("Number of Reactions : ", len(reaction))

        print("----------------------------- list of Calculation")
        reaction_str_len = max(len(r) for r in reaction)

        partial_xsec = []
        for i in range(num_cal):
            partial_sum_xsec = 0.0
            for j in range(len(data_matrix[i])):
                theta = np.radians(angle[j])
                d_theta = np.radians(angle_step)
                phi = 2 * np.pi
                partial_sum_xsec += data_matrix[i][j] * np.sin(theta) * d_theta * phi
            partial_xsec.append(partial_sum_xsec)

            pos = title[i].find(")")
            print(f"{reaction[i]:>{reaction_str_len + 3}} | {title[i][pos + 1:]} | Xsec({angle_min:3.0f}-{angle_max:3.0f} deg) : {partial_sum_xsec:.6f} mb")
        print("---------------------------------------------------")

        # ================================== save *.Ex.txt
        save_ex_name = f"{read_file[:-4]}.Ex.txt"
        print("Output : ", save_ex_name)
        with open(save_ex_name, 'w') as file_ex:
            file_ex.write("//generated_by_ExtractXSec.h____Ex____Xsec(4pi)____SF____sigma\n")
            for i in range(num_cal):
                file_ex.write(f"{Ex[i]:9.5f}     {partial_xsec[i]:9.5f}  1.0  0.000\n")
            file_ex.write("#=====End_of_File\n")

        # ================================== save file.Xsec.txt
        save_file_name = f"{read_file[:-4]}.Xsec.txt"
        print("Output : ", save_file_name)
        with open(save_file_name, 'w') as file_out:
            for r in reaction:
                file_out.write(f"#{r}\n")

            file_out.write("Angel\t")
            for t in title:
                file_out.write(f"{t:19}")
            file_out.write("\n")

            for i in range(len(angle)):
                file_out.write(f"{angle[i]:8.3f}\t")
                for j in range(num_cal):
                    file_out.write(f"{data_matrix[j][i]:19.6f}")
                file_out.write("\n")

    except Exception as e:
        print(f"An error occurred: {e}")

def is_float(value):
    try:
        float(value)
        return True
    except ValueError:
        return False