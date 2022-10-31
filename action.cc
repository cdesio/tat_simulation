#include "action.hh"

MyActionInitialization::MyActionInitialization()
{
}

MyActionInitialization::~MyActionInitialization()
{
}

void MyActionInitialization::Build() const
{
    MyPrimaryGenerator *generator = new MyPrimaryGenerator();
    SetUserAction(generator);

    // writing stuff to file
    MyRunAction *runAction = new MyRunAction();
    // using our own run action
    SetUserAction(runAction);

    // adding user defined event action
    MyEventAction *eventAction = new MyEventAction(runAction);
    SetUserAction(eventAction);

    // stepping action
    MySteppingAction *steppingAction = new MySteppingAction(eventAction);
    SetUserAction(steppingAction);
}