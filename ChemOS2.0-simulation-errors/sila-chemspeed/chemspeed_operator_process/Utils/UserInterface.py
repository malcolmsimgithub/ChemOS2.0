from .FileHandling import load_json, save_json
from .Timestamps import timestamp_time
import time
import telegram
import tkinter as tk
from tkinter import messagebox, scrolledtext
import socket



class UserInterface:
    def __init__(self, silasocket, telegram_comm=False, threads=("Main"), **kwargs):
        """Initializes a graphical user interface based on the tkinter library.
        Initializes the instance(s) of defined communicator(s), allowing for sending messages to a user.

        Currently implemented:
            - "print": all communication (output and input) via console output/input
            - "telegram": output communication via console and Telegram, input communication via console

        Parameters:
            telegram_comm (bool): if to send messages via Telegram
            threads (tuple or list): Collection of threads which are performed and logged individually
            kwargs (dict): communicator-specific kwargs

        Sets the following attributes:
            self.gui (GraphicalInterface)
            self.communicators (list): List of all active communicators.
        """
        if type(threads) is str:
            self.gui = GraphicalInterface(threads=[threads])
        else:
            self.gui = GraphicalInterface(threads=threads)
        self.communicators = [self.gui, PrintCommunicator(), Sila2Communicator(silasocket)]

        if telegram_comm:
            self.communicators.append(TelegramCommunicator(**kwargs))


    def activate_gui(self):
        """Activates the graphical user interface.
        ATTENTION: Blocks the main thread.

        Parameters:

        Returns:
            None
        """
        self.gui.run_main_loop()


    def communicate(self, message, thread="all", **kwargs):
        """Communicates a message via all active communicators.

        Parameters:
            message (str): Message to be communicated.
            thread (str): Name of the active thread

        Returns:
            None
        """
        for communicator in self.communicators:
            communicator.communicate(message, thread)


    def request_confirmation(self, message, thread="all", title="Confirmation Required", **kwargs):
        """Communicates a message via all active communicators. Requests confirmation in the Graphical User Interface.

        Parameters:
            message (str): Message to be communicated.
            title (str): Title of the popup window / prefix of the communicated message.
            thread (str): Name of the active thread

        Returns:
            None
        """
        time.sleep(5)
        self.communicate(f"{thread} - {title}: {message}", thread)
        self.gui.request_confirmation(message=message, title=title, thread=thread)



class PrintCommunicator:
    def __init__(self):
        """Initializes the print communicator (i.e. communicating all messages via console input / output)."""
        pass


    @staticmethod
    def communicate(message, thread):
        """Prints out the respective message to communicated, highlighted from regular console output.

        Parameters:
            message (str): Message to be communicated.
            thread (str): Name of the active thread

        Returns:
            None
        """
        print("")
        print("".rjust(len(message), "#"))
        print(f"{thread}: {message}")
        print("".rjust(len(message), "#"))


    def request_confirmation(self, message, thread):
        self.communicate(message, thread)
        while True:
            response = input("Please confirm (y/n)!")
            if response == "y":
                print("")
                return
            else:
                print("")

class Sila2Communicator:
    def __init__(self, silasocket):
        """Initializes the print communicator (i.e. communicating all messages via console input / output)."""

        self.silasocket= silasocket

    def communicate(self, message, thread):

        HEADER_LENGTH = 10

        # Prepare username and header and send them
        # We need to encode username to bytes, then count number of bytes and prepare header of fixed size, that we encode to bytes as well
        message = (thread+": " + message).encode('utf-8')
        message_header = f"{len(message):<{HEADER_LENGTH}}".encode('utf-8')
        self.silasocket.send(message_header+ message)


class GraphicalInterface:
    def __init__(self, threads):
        """Initializes a Graphical Interface Using the tkinter Library.
            - sends messages into a user interface "log window"
            - buttons for interrupting operation, unmounting the currently loaded tool, shutting down manager
            - can create popup messages to confirm
            - further functionalities could be included (e.g. via ttk, simpledialog, filedialog, ...)

        Parameters:
            threads (tuple or list): Collection of threads which should get individual columns

        ATTENTION: BLOCKS THE MAIN THREAD
        """
        self.root = tk.Tk()
        #self.root.state("zoomed")  # full active screen size
        self.root.title("ChemSpeed Operator")
        self.top = tk.Frame(self.root)
        self.main = tk.Frame(self.root)
        self.bottom = tk.Frame(self.root)
        self.threads = threads
        self.text = dict()
        self.confirmvar = 0
        for thread in threads:
            self.text[thread] = scrolledtext.ScrolledText(self.root)
        self.root.protocol("WM_DELETE_WINDOW", self.disable_close_button)  # substitutes the "close" [X] button
        self.pause_button = tk.Button(self.root, height=5, width=20, text="Pause Operation",
                                      command=(lambda: self.set_status("interrupt")))
        self.unmount_button = tk.Button(self.root, height=5, width=20, text="Unmount All",
                                        command=(lambda: self.set_status("emergency")))
        self.resume_button = tk.Button(self.root, height=5, width=20, text="Resume Operation",
                                       command=(lambda: self.set_status("restart")))
        self.shutdown_button = tk.Button(self.root, height=5, width=20, text="Shut Down Manager",
                                         command=(lambda: self.set_status("shutdown")))
        self.confirmButton = tk.Button(self.root, height=5, width=20, text="Confirm",
                                         command=(lambda: self.confirm()))



        self.status = None 
    def confirm(self):
        if self.status == "confirm required":
            self.confirmvar=1
        else:
            self.confirmvar=0

    def run_main_loop(self):
        """Function that opens the frame and runs through the main loop of the Tkinter object.
        IMPORTANT: Needs to be executed in the main thread!

        Parameters:

        Returns:
            None
        """
        self.bottom.pack(side="bottom")
        self.top.pack(side="top", fill="both", expand=True)
        for thread in self.threads:
            self.text[thread].pack(in_=self.top, expand=True, fill="both", side="left")
        self.pause_button.pack(in_=self.bottom, side="left")
        self.unmount_button.pack(in_=self.bottom, side="left")
        self.resume_button.pack(in_=self.bottom, side="left")
        self.shutdown_button.pack(in_=self.bottom, side="left")
        # self.confirmButton.pack(in_=self.bottom, side="left")
        self.root.mainloop()


    def communicate(self, message, thread):
        """Prints a message on the screen log.

        Parameters:
            message (str): Message to print
            thread (str): Which thread window to print to

        Returns:
            None
        """
        if thread == "all":
            thread = self.threads
        else:
            thread = [thread]

        for thr in thread:
            self.text[thr].insert("end", f"{timestamp_time()}   {message} \n")


    def request_confirmation(self, message, thread=None, title="Confirmation required!", return_reply=False):
        """Creates a popup window for the user to confirm the message.
        Freezes until the confirmation is given.

        Parameters:
            message (str): message to be printed
            thread (str or None): Name of the thread to be printed to
            title (str): Title of the popup window.
            return_reply (bool): Whether to return the user input as boolean or to wait until the input is True.

        Returns:
            None
        """
        self.status = "confirm required"
        if return_reply:
            return messagebox.askyesno(title, message)

        else:

            self.communicate(message,thread=thread)

            # while self.confirmvar == 0:
            #     time.sleep(1)
            
            # self.confirmvar = 0
            
            # self.communicate(f"Confirmed.", thread)
            self.status = "restart"
            
            # while True:
            #     confirmation = messagebox.askyesno(title, message)
            #     if confirmation:
            #         if thread:
            #             self.communicate(f"{message} (confirmed).", thread)
            #             return
            #         else:
            #             for thread in self.threads:
            #                 self.communicate(f"{message} (confirmed).", thread)
            #             return
            #     else:
            #         time.sleep(1)


    def set_status(self, status):
        """Sets the emergency status as True.

        Parameters:
            status (str or None): Target status

        Returns:
            None
        """
        if status in ("restart", "shutdown"):

            while True:

                # self.request_confirmation(f"Are you sure that the status should be set to {status}?", thread=None)
                # self.status = status
                confirmation = messagebox.askyesno("Confirmation required!", 
                f"Are you sure that the status should be set to {status}?")
                if confirmation:
                    self.status = status
                    self.communicate(f"Confirmed! setting status to {status}", thread='all')
                    return
                else:
                    return

        else:
            self.status = status
            return


    def disable_close_button(self):
        """Empty function that is just to be called when trying to hit the "close" button of the window.
        Prevents users from accidently closing the window.

        Parameters:

        Returns:
            None
        """
        pass


class TelegramCommunicator:
    def __init__(self, path_internal, recipient=None):
        """Sets up the Telegram bot to communicate via Telegram.

        Parameters:
            path_internal (pathlib.Path): Path where the "Utils" folder is located.
            recipient (str): Username of the active recipient. Must be registered in registered_users.json

        Sets the following attributes:
            self.path_internal (pathlib.Path): Path to the internal "Telegram_Data" folder
            self.token (str): Identifier token of the Telegram bot.
            self.bot (telegram.Bot): Telegram Bot
            self.updater (telegram.ext.Updater): Handler for incoming messages
            self.registered_users (dict): Dictionary of registered users and their chat IDs {username: chat_id}
            self.active_chats (dict): Dictionary of active chats {chat_id: {"username": username, "status": status}}
            self.recipient (str): Chat ID of the current standard recipient of communications.
        """
        self.path_internal = path_internal / "Utils" / "Telegram_Data"

        self.token = "5077799396:AAGEMZmRKR3_D6qLm8P4Vw50LY6o3vWZByk"
        self.bot = telegram.Bot(self.token)
        self.updater = telegram.ext.Updater(token=self.token, use_context=True)
        self.updater.dispatcher.add_handler(telegram.ext.MessageHandler(telegram.ext.filters.Filters.update, callback=self.incoming_message))
        self.updater.start_polling()

        self.registered_users, self.active_chats = self.load_existing_users()

        if recipient:
            try:
                self.recipient = self.registered_users[recipient]
                self.active_chats[self.recipient]["status"] = "active"
            except KeyError:
                print(f"The user name {recipient} has not been registered yet.")
                self.recipient = None
        else:
            self.recipient = None

    #############################################################
    # Low-level methods to be called by the UserInterface class #
    #############################################################

    def communicate(self, message):
        """Checks if recipient is specified, and sends the message.

        Parameters:
            message (str): Message to be sent.

        Returns:
            None
        """
        if self.recipient:
            self.send_message(message, self.recipient)

    def request_confirmation(self, message):
        """Currently: Confirmations only possible via command line confirmation. Only the notification message is sent.

        Parameters:
            message (str): Message to be sent.

        Returns:
            None
        """
        self.communicate(message)

    ###################################################################
    # Low-level methods to manage incoming and outgoing communication #
    ###################################################################

    def send_message(self, message, recipient, **kwargs):
        """Sends a given message to a recipient (identified by their chat ID).

        Parameters:
            message (str): Message to be sent
            recipient (int): Chat ID of the recipient.
            kwargs (dict): Further keyword arguments

        Returns:
            None
        """
        self.bot.sendMessage(text=message, chat_id=recipient, **kwargs)

    def incoming_message(self, update, context):
        """Handler for incoming messages (as coming from updater/dispatcher – ext. message handler).
        Calls specific methods depending on the sender of the message.

        Parameters:
            update (telegram.Update): Update sent by the telegram.updater
            context (telegram.Context): Context of the update

        Returns:
            None
        """
        chat_id = update.effective_chat.id
        sender_status = self.check_sender_status(chat_id)

        if sender_status == "unknown":
            if update.message.text == "REGISTER":
                self.register_new_user(chat_id)

        elif sender_status == "incomplete_pw":
            self.check_password(chat_id, update)

        elif sender_status == "incomplete_un":
            self.complete_registration(chat_id, update)

        elif sender_status == "registered":
            self.message_from_user(chat_id, update)

        elif sender_status == "active":
            self.message_from_active_user(chat_id, update)

    ##############################################################
    # High-level internal methods to process specific threads #
    ##############################################################

    def register_new_user(self, chat_id):
        """Registers a new user as incomplete in self.active_chats and requests password.

        Parameters:
            chat_id (int): Chat ID of the new user

        Returns:
            None
        """
        self.send_message("Welcome as a new potential ChemSpeed User!", chat_id)
        self.send_message("To register, please enter the password.", chat_id)
        self.active_chats[chat_id] = {"username": None, "status": "incomplete_pw"}

    def check_password(self, chat_id, update):
        """Checks password given by pre-registered user and asks for user name (if correct).

        Parameters:
            chat_id (int): Chat ID of the new user
            update (telegram.Update): Incoming message / update, as provided by telegram.updater
        """
        if update.message.text == "chemspeed_aspuru":
            self.send_message("The password is correct.", chat_id)
            self.send_message("Please choose your internal username.", chat_id)
            self.active_chats[chat_id]["status"] = "incomplete_un"
        else:
            self.send_message("Wrong passphrase.", chat_id)

    def complete_registration(self, chat_id, update):
        """Completes registration of the user (registration in self.registered_users and the corresponding file). Sets the user name.

        Parameters:
            chat_id (int): Chat ID of the new user
            update (telegram.Update): Incoming message / update, as provided by telegram.updater
        """
        self.active_chats[chat_id]["username"] = update.message.text
        self.active_chats[chat_id]["status"] = "registered"
        self.send_message(f"Success! The user {update.message.text} was registered.", chat_id)
        self.registered_users[update.message.text] = chat_id
        self.save_existing_users()

    def load_existing_users(self):
        """Loads all previously registered users from the internal registered_users.json file.

        Returns:
            user_dict (dict): Dictionary of all users {username: chat_id}
            active_chats (dict): Dictionary of all registered user chats and their statuses {chat_id: {"username": username, "status": "registered"}}
        """
        user_dict = load_json(self.path_internal / "registered_users.json")
        active_chats = dict()

        for user in user_dict:
            active_chats[user_dict[user]] = {"username": user, "status": "registered"}

        return user_dict, active_chats

    def save_existing_users(self):
        """Saves all currently registered users to the internal registered_users.json file.
        Can be used for updating the user register.

        Returns:
            None
        """
        save_json(self.registered_users, self.path_internal / "registered_users.json")

    def message_from_user(self, user_id, update):
        """Handles incoming messages from registered users (currently: print the message).

        Parameters:
            user_id (int): Chat ID of the user
            update (telegram.Update): Incoming message / update, as provided by telegram.updater
        """
        print(f"{self.active_chats[user_id]['username']} sent the following message:")
        print(update.message.text)
        print("")

    def message_from_active_user(self, user_id, update):
        """Handles incoming messages from the currently active user (currently: no specific implementation for active user).

        Parameters:
            user_id (int): Chat ID of the user
            update (telegram.Update): Incoming message / update, as provided by telegram.updater
        """
        self.message_from_user(user_id, update)

    def check_sender_status(self, user_id):
        """Checks the user status for the user of an incoming message.

        Parameters:
            user_id (int): Chat ID of the user

        Returns:
            status (str): User status ("unknown", "incomplete_pw", "incomplete_un", "registered", "active")
        """
        if user_id not in self.active_chats:
            return "unknown"
        else:
            return self.active_chats[user_id]["status"]
