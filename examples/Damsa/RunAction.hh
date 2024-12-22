#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "SteppingActionDMG4.hh"


class RunAction : public G4UserRunAction
{
    public:
        RunAction();
        virtual ~RunAction();

        //virtual void BeginOfRunAction(const G4Run* run);
        virtual void EndOfRunAction(const G4Run* run);
        virtual void BeginOfRunAction(const G4Run* run);
        
    private:
        SteppingActionDMG4* fSteppingAction;
        // stepping action.hh에서 변수 설정한 걸로 steppingAction.cc에서 데이터를 쌓고 여기서 stepAction.hh의 객체 타입을 불러옴. 
        // .cc 파일에서 한번만 사용할거면 거기서 한번만 설정하면 되는데 그게 아니라 계속 사용할거라면 
        // 그냥 여기서 전역변수 하나 만들어 놓고 이거 맨 처음에 초기화 시켜놓고 계속 쓰는거라 함 
};

#endif