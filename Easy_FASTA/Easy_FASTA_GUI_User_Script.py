# Libraries
import re
import PySimpleGUI as sg
import os
import pyautogui


# Define a function containing all the code structure. Needed for creating the button menu
def fasta_search():

    # Select a theme
    sg.theme("DarkAmber")

    # Set the menu to show up upon clicking on the "Search by" button
    but_menu = ["Unused", ["Gene tag", "Locus tag", "Main tag", "Protein ID", "Protein name", "Niben261", "NbLAB360"]]

    # List of keys to clear after the end of the process
    clear_keys = ["-TEXT_ID-", "-TEXT_FASTA-", "-TEXT_FOLDER-", "-SAVE-", "-ONELINE-"]

    # Create a window layout
    layout = [[sg.Image(filename="Easy_FASTA_Logo.png")],

              [sg.Text(text="IDs file:", tooltip="A \".txt\" file with IDs", pad=((25, 0), (10, 0)))],
              [sg.In(size=(29, 1), pad=((25, 0), (0, 0)), justification="c", key="-TEXT_ID-"),
               sg.FileBrowse(enable_events=True, key="-ID-", size=(8, 1), pad=((20, 0), (0, 0)),
                             file_types=(("Text Files", "*.txt"),))],

              [sg.Text(text="FASTA file:", pad=((25, 0), (0, 0)))],
              [sg.In(size=(29, 1), pad=((25, 0), (0, 0)), justification="c", key="-TEXT_FASTA-"),
               sg.FileBrowse(enable_events=True, key="-FASTA-", size=(8, 1), pad=((20, 0), (0, 0)),
                             file_types=(("*.txt", "*.txt"), ("*.fna", "*.fna"), ("*.faa", "*.faa"), ("*.fasta", "*.fasta")))],

              [sg.Text(text="Save directory:", tooltip="Where to save the output file", pad=((25, 0), (0, 0)))],
              [sg.In(size=(29, 1), pad=((25, 0), (0, 0)), justification="c", key="-TEXT_FOLDER-"),
               sg.FolderBrowse(enable_events=True, key="-FOLDER-", size=(8, 1), pad=((20, 0), (0, 0)))],

              [sg.Text(text="Output file:", tooltip="The name of the output file", pad=((25, 0), (0, 0))),
               sg.Checkbox(text="One line sequence?", tooltip="Whether to save each target sequence in a single line",
                           enable_events=True, key="-ONELINE-")],
              [sg.In(size=(29, 1), pad=((25, 0), (0, 0)), justification="c", enable_events=True, key="-SAVE-")],

              [sg.ButtonMenu(button_text="Search by", menu_def=but_menu, key="-SEARCH-",
                             size=(8, 1), pad=((55, 0), (30, 0)), tooltip="The type of ID in the \"IDs file\"."),
               sg.Button(button_text="Run", size=(8, 1), pad=((20, 0), (30, 0))),
               sg.Button(button_text="Clear", size=(8, 1), pad=((20, 0), (30, 0)))]]

    # Create a window
    window = sg.Window("Easy FASTA", layout, font="Arial 11",
                       icon="Easy_FASTA_Icon.ico", size=(450, 450),
                       resizable=True)
                       

    # Open the window
    while True:
        event, values = window.read()
        if event == sg.WIN_CLOSED:
            break

        # Open Input ID file
        ID_count = 0
        IDs = None
        IDs_file = values["-ID-"]
        if ".txt" not in IDs_file:
            ID_count = ID_count + 1
        if len(IDs_file) > 0:
            IDs = open(IDs_file)

        # Open FASTA file
        FASTA = None
        FASTA_II = None
        FASTA_count = 0
        FASTA_file = values["-FASTA-"]
        extensions = [".txt", ".fna", ".faa", ".fasta"]  # Only these extensions will be accepted
        for element in extensions:
            if element not in FASTA_file:
                FASTA_count = FASTA_count + 1
        if len(FASTA_file) > 0:
            FASTA = open(FASTA_file)  # It will only be used for obtaining the ID
            FASTA_II = open(FASTA_file)  # It will be used for obtaining the sequences

        # Output file name
        fsave = values["-SAVE-"]
        if ".txt" not in fsave and len(fsave) > 0:
            fsave = fsave + ".txt"
        empty = re.findall(".", fsave)

        # Save directory
        user_path = values["-FOLDER-"]
        fsave = os.path.join(user_path, fsave)
        path = re.findall(".:", user_path)  # Only accept actual directories. "path" is a list.

        # Search by method
        search = values["-SEARCH-"]

        # Clear the input fields
        if event == "Clear":
            for key in clear_keys:
                window[key]("")

        # Error Messages
        if event == "Run" and len(values["-TEXT_ID-"]) == 0:
            sg.PopupError("No IDs file was selected.", title="ERROR MESSAGE")
            continue
        elif event == "Run" and ID_count == 1:
            sg.PopupError("The IDs file must contain the \".txt\" extension.", title="ERROR MESSAGE")
            continue
        elif event == "Run" and len(values["-TEXT_FASTA-"]) == 0:
            sg.PopupError("No FASTA file was selected.", title="ERROR MESSAGE")
            continue
        elif event == "Run" and FASTA_count == 4:
            sg.PopupError("The FASTA file must contain one of the following extensions:",
                          "                          ""\".txt\", \".faa\", \".fna\" or \".fasta\".", title="ERROR MESSAGE")
            continue
        elif event == "Run" and len(values["-TEXT_FOLDER-"]) == 0:
            sg.PopupError("No save directory was selected.", title="ERROR MESSAGE")
            continue
        elif event == "Run" and len(path) == 0:
            sg.PopupError("Invalid save directory.", title="ERROR MESSAGE")
            continue
        elif event == "Run" and len(empty) == 0:
            sg.PopupError("Output file name not specified.", title="ERROR MESSAGE")
            continue
        elif event == "Run" and search is None:
            sg.PopupError("No \"search by\" method was specified.", title="ERROR MESSAGE")
            continue

        # Append the IDs from the input file to a list
        if event == "Run":
            # Determine the prefixes (especially useful for the "Protein name" and "Gene tag" options, as they can be
            # similar to any other word/part of word in the line, thus not making the program to retrieve them properly
            prefix = ""
            if values["-SEARCH-"] == "Main tag" or values["-SEARCH-"] == "Niben261" or values["-SEARCH-"] == "NbLAB360":
                prefix = ""
            if values["-SEARCH-"] == "Protein ID":
                prefix = "protein_id="
            if values["-SEARCH-"] == "Protein name":
                prefix = "protein="
            if values["-SEARCH-"] == "Locus tag":
                prefix = "locus_tag="
            if values["-SEARCH-"] == "Gene tag":
                prefix = "gene="

            # Append IDs into a list
            lst_IDs = list()
            lst_IDs_II = list()
            lst_IDs_fixed = list()
            lst_tag_check = list()

            for lines in IDs:
                lines = lines.strip()

                if values["-SEARCH-"] == "Niben261" or values["-SEARCH-"] == "NbLAB360":
                    tag = re.findall("(\S+\.\S)", lines)
                    if len(tag) > 0:
                        lst_IDs.append(prefix + tag[0])
                        lst_IDs_II.append(prefix + tag[0])
                        lst_IDs_fixed.append(tag[0])
                        lst_tag_check.append(tag[0])

                else:        
                    lst_IDs.append(prefix + lines)
                    lst_IDs_II.append(prefix + lines)
                    lst_IDs_fixed.append(lines)
                    lst_tag_check.append(lines)

            lst_IDs = list(dict.fromkeys(lst_IDs))  # Remove any duplicate values to avoid duplicate outputs
            lst_IDs = [x.lower() for x in lst_IDs] # Convert input IDs to lowercase
            lst_IDs_II = list(dict.fromkeys(lst_IDs_II))  # Remove any duplicate values to avoid duplicate outputs
            lst_IDs_II = [x.lower() for x in lst_IDs_II] # Convert input IDs to lowercase
            duplicates = len(lst_IDs_fixed) - len(lst_IDs)

            # Append the lines with IDs from the FASTA file to a list
            lst_FASTA = list()
            for lines in FASTA:
                lines = lines.strip()
                if not lines.startswith(">"):
                    continue
                lst_FASTA.append(lines)

            lst_FASTA = [x.lower() for x in lst_FASTA] # Convert input IDs to lowercase

            # Sort the user desired IDs based on "lst_FASTA" to make the parsing process much faster
            lst_IDs_sorted = list()
            tag = None
            for element in lst_FASTA:

                if values["-SEARCH-"] == "Main tag":
                    tag = re.findall("^>lcl\|(\S+)", element)
                    if len(tag) == 0:
                        continue

                if values["-SEARCH-"] == "Protein ID":
                    tag = re.findall("(protein_id=\S+)\]", element)
                    if len(tag) == 0:
                        continue

                if values["-SEARCH-"] == "Protein name":
                    tag = re.findall("(protein=.+?)\]", element)
                    if len(tag) == 0:
                        continue

                if values["-SEARCH-"] == "Locus tag":
                    tag = re.findall("(locus_tag=\S+)\]", element)
                    if len(tag) == 0:
                        continue

                if values["-SEARCH-"] == "Gene tag":
                    tag = re.findall("(gene=.+?)\]", element)
                    if len(tag) == 0:
                        continue

                if values["-SEARCH-"] == "Niben261" or values["-SEARCH-"] == "NbLAB360":
                    tag = re.findall("^>(\S+\.\S)\s", element)
                    if len(tag) == 0:
                        continue
                
                
                for element_II in lst_IDs:
                    if element_II == tag[0]:  # Substitutes the "re.fullmatch()" method. Only an exact match will be accepted
                        lst_IDs_sorted.append(element_II)

            # Parse the FASTA file to get the desired sequences
            count = 0
            ID_lst = list()
            seq_lst = list()
            seq_lst_fixed = list()
            for lines in FASTA_II:
                lines = lines.strip()
                if lines.startswith(">") and count < len(lst_IDs_sorted):
                    lines = lines.lower()
                    if values["-SEARCH-"] == "Main tag":
                        tag = re.findall("^>lcl\|(\S+)", lines)
                    if values["-SEARCH-"] == "Protein ID":
                        tag = re.findall("(protein_id=\S+)\]", lines)
                    if values["-SEARCH-"] == "Protein name":
                        tag = re.findall("(protein=.+?)\]", lines)
                    if values["-SEARCH-"] == "Locus tag":
                        tag = re.findall("(locus_tag=\S+)\]", lines)
                    if values["-SEARCH-"] == "Gene tag":
                        tag = re.findall("(gene=.+?)\]", lines)
                    if values["-SEARCH-"] == "Niben261" or values["-SEARCH-"] == "NbLAB360":
                        tag = re.findall("^>(\S+\.\S)\s", lines)

                    if len(tag) == 0:
                        ID_lst.clear()  # This won't allow getting sequences that do not contain the type of tag
                        continue
                    if lst_IDs_sorted[count] != tag[0]:  # A much better way than using "in"
                        ID_lst.clear()  # This won't allow for getting sequences that do not contain the type of tag
                        continue
                    else:
                        ID_lst.append(lines)
                        if len(ID_lst) > 1:
                            ID_lst.remove(
                                ID_lst[0])  # Remove only one element so that len(ID_lst) is always 1 if there's
                            # a match
                        count += 1
                        if count == len(lst_IDs_sorted):
                            # Fix the count value when the the app reaches the last sequence. Useful for the
                            # progress bar.
                            count = len(lst_IDs_sorted) - 1

                # If there's a match, get the sequence and put it into a list -> "seq_lst"
                if len(ID_lst) == 1:
                    seq_lst.append(lines)
                    seq_lst_fixed.append(lines)
                    if lines.startswith(">"):  # Remove found IDs below here rather than the commented above
                        for element_III in lst_IDs_II:
                            if element_III in lines:
                                lst_IDs_II.remove(element_III)  # Keep the not found IDs to show later
                                if not sg.one_line_progress_meter("Progress", current_value=count,
                                                                  max_value=len(lst_IDs_sorted)):
                                    break

            # Save the sequences into a file
            save_count = 0
            lst_ID = list()
            lst_seq = list()
            count_seq = 0
            for element in seq_lst:
                if values["-ONELINE-"] is True:
                    if element.startswith(">"):
                        lst_ID.append(element)

                    if len(lst_IDs_sorted) == 1:
                        with open(fsave, "a") as f:
                            print(seq_lst[0], file=f)  # The line with IDs
                            print(seq_lst[1:], sep="", file=f)  # The sequence
                            save_count = len(seq_lst_fixed)

                    if count_seq == len(lst_IDs_sorted) - 1:
                        seq_lst = seq_lst[save_count:]  # The start position of the last sequence is equals to
                        # save count
                        with open(fsave, "a") as f:
                            print(lst_ID[0], file=f)
                            print(*seq_lst, sep="", file=f)
                            save_count = len(seq_lst_fixed)

                    else:
                        if len(lst_ID) > 1:
                            with open(fsave, "a") as f:
                                print(lst_ID[0], file=f)
                                print(*lst_seq, sep="", file=f)
                            lst_seq.clear()
                            lst_ID.remove(lst_ID[0])
                            count_seq += 1

                        if not element.startswith(">"):
                            lst_seq.append(element)
                        save_count += 1

                else:
                    with open(fsave, "a") as f:
                        print(element, file=f)
                        save_count += 1

                # Set a progress save bar
                if not sg.one_line_progress_meter("Progress", current_value=save_count,
                                                  max_value=len(seq_lst_fixed)):
                    break  # The break will allow for the user to use the "Cancel" button

            # End of the process message
            if count == len(lst_IDs_sorted) - 1:
                sg.PopupOK("Number of input IDs (nID): " + str(len(lst_IDs_fixed)),
                           "Number of retrieved sequences (nSeq): " + str(count + 1),
                           "Number of duplicates: " + str(duplicates),
                           "Tip: If nSeq > nID, at least one input ID repeats in the FASTA file",
                           title="End of the process")

            else:  # In case some IDs were not found
                sg.PopupOK("Number of input IDs (nID): " + str(len(lst_IDs_fixed)),
                           "Number of retrieved sequences (nSeq): " + str(count),
                           "Number of duplicates: " + str(duplicates),
                           "Tip: If nSeq > nID, at least one input ID repeats in the FASTA file",
                           title="End of the process")

            # PopUp for not found IDs
            if len(lst_IDs_II) > 0:
                not_found = None
                if len(lst_IDs_II) == 1:
                    not_found = sg.PopupYesNo(str(len(lst_IDs_II)) + " ID was not found. Do you wish to save it?",
                                              title="No match ID")
                if len(lst_IDs_II) > 1:
                    not_found = sg.PopupYesNo(str(len(lst_IDs_II)) + " IDs were not found. Do you wish to save them?",
                                              title="No match IDs")
                if not_found == "Yes":
                    while True:
                        save_name = sg.PopupGetText("Save file name:", title="Error file")
                        if save_name is None:
                            break
                        if len(save_name) == 0:
                            continue
                        else:
                            if ".txt" not in save_name:
                                save_name = save_name + ".txt"
                            while True:
                                save_path = sg.PopupGetFolder("Choose the save directory:", title="Save directory")
                                if save_path is None or len(save_path) > 0:
                                    break
                                else:
                                    continue
                            if save_path is None:
                                break
                            else:
                                save_name = os.path.join(save_path, save_name)
                                for element in lst_IDs_II:
                                    with open(save_name, "a") as f:
                                        print(element, file=f)

                                # As the user whether to access the save directory of the error file
                                if "/" in save_path:
                                    save_path = save_path.replace("/", "\\")
                                cd = sg.PopupYesNo("Do you want to access the output error file directory?",
                                                   title="Access save directory")
                                if cd == "Yes":
                                    pyautogui.press("winleft")
                                    pyautogui.write(save_path)
                                    pyautogui.press("enter")
                                    pyautogui.press("enter")
                                break

            # Ask the user whether to access the save directory of the output sequence file
            if len(seq_lst_fixed) > 0:
                if "/" in user_path:
                    user_path = user_path.replace("/", "\\")
                cd = sg.PopupYesNo("Do you want to access the output FASTA file directory?",
                                   title="Access save directory")
                if cd == "Yes":
                    pyautogui.press("winleft")
                    pyautogui.write(user_path)
                    pyautogui.press("enter")
                    pyautogui.press("enter")

    window.close()

fasta_search()
