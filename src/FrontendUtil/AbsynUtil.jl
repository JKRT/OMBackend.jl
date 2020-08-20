module AbsynUtil

import Absyn
using MetaModelica
using Absyn
const dummyParts = PARTS(nil, nil, nil, nil, NONE())::ClassDef
const dummyInfo = SOURCEINFO("", false, 0, 0, 0, 0, 0.0)::Info
const dummyProgram = PROGRAM(nil, TOP())::Program

end #=End AbsynUtil=#
